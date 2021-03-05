/**
MIT License

Copyright (c) 2017 DECTRIS Ltd.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include <dectris/neggia/data/H5Superblock.h>
#include <dectris/neggia/data/JenkinsLookup3Checksum.h>
#include <dectris/neggia/user/Dataset.h>
#include <dectris/neggia/user/H5File.h>
#include <gtest/gtest.h>
#include <numeric>
#include <type_traits>

struct DatasetValues {
    std::string entry;
    size_t frame;
    std::vector<size_t> dim;
    uint64_t valuesum;
    uint32_t checksum;
};

template <typename FloatParameter,
          typename IntegerParameter,
          typename PixelType>
struct ExpectedValues {
    uint8_t superblock_version;
    IntegerParameter width;
    IntegerParameter height;
    FloatParameter x_pixel_size;
    FloatParameter y_pixel_size;
    uint64_t pixel_mask_valuesum;
    uint32_t pixel_mask_checksum;
    std::vector<DatasetValues> datasets;
};

template <typename ValueType>
void CheckHdf5SingleValue(const std::string& filename,
                          const std::string& h5path,
                          const ValueType& expected) {
    Dataset ds(H5File(filename), h5path);
    ASSERT_TRUE(ds.dim().empty());
    ASSERT_FALSE(ds.isChunked());
    if (std::is_same<ValueType, float>::value ||
        std::is_same<ValueType, double>::value)
    {
        ASSERT_EQ(ds.dataTypeId(), 1);
    } else {
        ASSERT_EQ(ds.dataTypeId(), 0);
    }
    ASSERT_EQ(ds.dataSize(), sizeof(ValueType));
    ValueType val;
    ds.read(&val);
    ASSERT_EQ(val, expected);
};

template <typename ValueType>
void CheckIntegerDataset(const std::string& filename,
                         const std::string& h5path,
                         size_t frame,
                         const std::vector<size_t>& dim,
                         uint64_t valuesum,
                         uint32_t checksum) {
    Dataset ds(H5File(filename), h5path);
    ASSERT_EQ(ds.dim(), dim);
    ASSERT_EQ(ds.dataSize(), sizeof(ValueType));
    size_t pixel_count;
    std::vector<size_t> chunkOffset;
    switch (ds.dim().size()) {
        case 2:
            ASSERT_EQ(frame, 0);
            pixel_count = ds.dim().at(0) * ds.dim().at(1);
            break;
        case 3:
            pixel_count = ds.dim().at(1) * ds.dim().at(2);
            chunkOffset = std::vector<size_t>({frame, 0, 0});
            break;
        default:
            throw std::runtime_error("dataset with dimensionality of " +
                                     std::to_string(ds.dim().size()) +
                                     " not supported.");
    }

    ValueType datasetArray[pixel_count];
    ds.read(datasetArray, chunkOffset);

    uint64_t sum = std::accumulate(datasetArray, datasetArray + pixel_count,
                                   uint64_t(0));
    ASSERT_EQ(sum, valuesum);

    uint32_t checkSumCalculated = JenkinsLookup3Checksum(std::string(
            (const char*)datasetArray, pixel_count * sizeof(ValueType)));
    ASSERT_EQ(checkSumCalculated, checksum);
}

template <typename FloatParameter,
          typename IntegerParameter,
          typename PixelType>
void CheckHdf5(
        const std::string& filename,
        const ExpectedValues<FloatParameter, IntegerParameter, PixelType>&
                expected) {
    H5File h5File(filename);
    H5Superblock superblock(h5File.fileAddress());
    ASSERT_EQ(superblock.version(), expected.superblock_version);
    CheckHdf5SingleValue<>(
            filename,
            "/entry/instrument/detector/detectorSpecific/x_pixels_in_detector",
            expected.width);
    CheckHdf5SingleValue<>(
            filename,
            "/entry/instrument/detector/detectorSpecific/y_pixels_in_detector",
            expected.height);
    CheckHdf5SingleValue<>(filename, "/entry/instrument/detector/x_pixel_size",
                           expected.x_pixel_size);
    CheckHdf5SingleValue<>(filename, "/entry/instrument/detector/y_pixel_size",
                           expected.y_pixel_size);
    CheckIntegerDataset<uint32_t>(
            filename, "/entry/instrument/detector/detectorSpecific/pixel_mask",
            0, {expected.height, expected.width}, expected.pixel_mask_valuesum,
            expected.pixel_mask_checksum);
    for (auto dataset : expected.datasets) {
        CheckIntegerDataset<PixelType>(filename, "/entry/data/" + dataset.entry,
                                       dataset.frame, dataset.dim,
                                       dataset.valuesum, dataset.checksum);
    }
};

TEST(TestDatasetEiger1, MasterOnlyBSLZ4) {
    CheckHdf5<float, uint32_t>(
            "h5-testfiles/dataset_eiger1_001/"
            "eiger1_testmode10_0datafiles_4images_bslz4_master.h5",
            ExpectedValues<float, uint32_t, uint32_t>{
                    0,
                    1030,
                    1065,
                    7.5e-5,
                    7.5e-5,
                    38372,
                    2854193483,
                    {},
            });
}

TEST(TestDatasetEiger2, BSLZ4) {
    CheckHdf5<double, uint64_t>(
            "h5-testfiles/dataset_eiger2_001/eiger2_master.h5",
            ExpectedValues<double, uint64_t, uint16_t>{
                    2,
                    1030,
                    1064,
                    7.5e-5,
                    7.5e-5,
                    47344,
                    3591651806,
                    {{"data_000001",
                      0,
                      {3, 1064, 1030},
                      71821117200,
                      1187733511}},
            });
}