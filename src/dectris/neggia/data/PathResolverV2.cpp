#include "PathResolverV2.h"
#include <assert.h>
#include <memory>
#include "H5BTreeVersion2.h"
#include "H5FractalHeap.h"
#include "H5LocalHeap.h"
#include "constants.h"

#include <iostream>

PathResolverV2::PathResolverV2(const H5ObjectHeader& root) : _root(root) {}

ResolvedPath PathResolverV2::resolve(const H5Path& path) {
    return resolvePathInHeader(_root, path);
}

ResolvedPath PathResolverV2::resolvePathInHeader(const H5ObjectHeader& in,
                                                 const H5Path& path) {
    H5ObjectHeader parentEntry(path.isAbsolute() ? _root : in);
    for (auto itemIterator = path.begin(); itemIterator != path.end();
         ++itemIterator)
    {
        auto item = *itemIterator;
        try {
            auto resolvedPath = findPathInObjectHeader(
                    parentEntry, item, H5Path(path, itemIterator + 1));
            if (resolvedPath.externalFile.get() != nullptr) {
                // the object is in a different file which must be opened by
                // upper layer
                return resolvedPath;
            }
            if ((itemIterator + 1) == path.end()) {
                return resolvedPath;
            }
            parentEntry = resolvedPath.objectHeader;
        } catch (const std::out_of_range&) {
            continue;
        }
    }
    auto output = ResolvedPath{};
    output.objectHeader = parentEntry;
    return output;
}

ResolvedPath PathResolverV2::findPathInLinkMsg(
        const H5ObjectHeader& parentEntry,
        const H5LinkMsg& linkMsg,
        const H5Path& remainingPath) {
    switch (linkMsg.linkType()) {
        case H5LinkMsg::SOFT: {
            H5Path targetPath(linkMsg.targetPath());
            return resolvePathInHeader(parentEntry, targetPath + remainingPath);
        } break;
        case H5LinkMsg::HARD: {
            H5Path targetPath(linkMsg.targetPath());
            return resolvePathInHeader(linkMsg.hardLinkObjectHeader(),
                                       targetPath + remainingPath);
        } break;
        case H5LinkMsg::EXTERNAL: {
            std::string targetFile = linkMsg.targetFile();
            H5Path targetPath(linkMsg.targetPath());
            auto output = ResolvedPath{};
            output.externalFile.reset(new ResolvedPath::ExternalFile{
                    targetFile, targetPath + remainingPath});
            return output;
        } break;
        default: {
            throw std::runtime_error("unknown link type" + linkMsg.linkType());
        }
    }
}

uint32_t PathResolverV2::getFractalHeapOffset(
        const H5LinkInfoMsg& linkInfoMsg,
        const std::string& pathItem) const {
    size_t btreeAddress = linkInfoMsg.getBTreeAddress();
    if (btreeAddress == H5_INVALID_ADDRESS) {
        throw std::out_of_range("Invalid address");
    }
    H5BTreeVersion2 btree(_root.fileAddress(), btreeAddress);
    H5Object heapRecord(_root.fileAddress(), btree.getRecordAddress(pathItem));
    return heapRecord.read_u32(5);
}

ResolvedPath PathResolverV2::findPathInObjectHeader(
        const H5ObjectHeader& parentEntry,
        const std::string pathItem,
        const H5Path& remainingPath) {
    for (size_t i = 0; i < parentEntry.numberOfMessages(); ++i) {
        auto msg = parentEntry.headerMessage(i);
        switch (msg.type) {
            case H5LinkMsg::TYPE_ID: {
                H5LinkMsg linkMsg(msg.object);
                if (linkMsg.linkName() != pathItem)
                    continue;
                return findPathInLinkMsg(parentEntry, linkMsg, remainingPath);
            } break;
            case H5LinkInfoMsg::TYPE_ID: {
                uint32_t heapOffset;
                H5LinkInfoMsg linkInfoMsg(msg.object);
                try {
                    heapOffset = getFractalHeapOffset(linkInfoMsg, pathItem);
                } catch (const std::out_of_range&) {
                    continue;
                }
                H5FractalHeap fractalHeap(_root.fileAddress(),
                                          linkInfoMsg.getFractalHeapAddress());
                H5LinkMsg linkMsg(fractalHeap.getHeapObject(heapOffset));
                assert(linkMsg.linkName() == pathItem);
                auto output = ResolvedPath{};
                output.objectHeader = linkMsg.hardLinkObjectHeader();
                return output;
            } break;
            default: {
                continue;
            }
        }
    }
    throw std::out_of_range("Not found");
}