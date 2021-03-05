#!/usr/bin/env python3
from dataclasses import dataclass
import hdf5plugin # must be imported before h5py
import h5py
import os


@dataclass
class DatasetDescriptor:
    master_file: str
    superblock_version: int
    data_entries: list
    pm_checksum: int
    pm_shape: tuple
    ff_checksum: int
    ff_shape: tuple
    data_checksum: int
    data_shape: tuple



def describe(filename):
    with open(filename, "rb") as f:
        header = f.read(8)
        version = f.read(1)
    if header != b"\211HDF\r\n\032\n":
        print(f"BAD HEADER")
    superblock_version = int(version[0])
    master = h5py.File(filename, "r")
    data_entries = [str(key) for key in master["entry"]["data"].keys()]
    if "pixel_mask" in master["/entry/instrument/detector/detectorSpecific"].keys():
        pm_checksum = lookup3(master["/entry/instrument/detector/detectorSpecific/pixel_mask"][()].tobytes())
        pm_shape = master["/entry/instrument/detector/detectorSpecific/pixel_mask"][()].shape
    else:
        pm_checksum = -1
        pm_shape = ()
    if "flatfield" in master["/entry/instrument/detector/detectorSpecific"].keys():
        ff_checksum = lookup3(master["/entry/instrument/detector/detectorSpecific/flatfield"][()].tobytes())
        ff_shape = master["/entry/instrument/detector/detectorSpecific/flatfield"][()].shape
    else:
        ff_checksum = -1
        ff_shape = ()
    data_checksum = lookup3(master[f"/entry/data/{data_entries[0]}"][()].tobytes())
    data_shape = master[f"/entry/data/{data_entries[0]}"][()].shape
    for i in range(1, len(data_entries)):
        data_chk = lookup3(master[f"/entry/data/{data_entries[i]}"][()].tobytes())
        if data_chk != data_checksum:
            print(f"ERROR: hdf5 with changing data entries (checksum)")
        data_s = master[f"/entry/data/{data_entries[i]}"][()].shape
        if data_s != data_shape:
            print(f"ERROR: hdf5 with changing data entries (shape)")
    print(str(DatasetDescriptor(
        master_file=filename,
        superblock_version=superblock_version,
        data_entries=data_entries,
        pm_checksum=pm_checksum,
        pm_shape=pm_shape,
        ff_checksum=ff_checksum,
        ff_shape=ff_shape,
        data_checksum=data_checksum,
        data_shape=data_shape,
    )))

def rot(x,k):
    return (((x)<<(k)) | ((x)>>(32-(k))))

def mix(a, b, c):
    a &= 0xffffffff; b &= 0xffffffff; c &= 0xffffffff
    a -= c; a &= 0xffffffff; a ^= rot(c,4);  a &= 0xffffffff; c += b; c &= 0xffffffff
    b -= a; b &= 0xffffffff; b ^= rot(a,6);  b &= 0xffffffff; a += c; a &= 0xffffffff
    c -= b; c &= 0xffffffff; c ^= rot(b,8);  c &= 0xffffffff; b += a; b &= 0xffffffff
    a -= c; a &= 0xffffffff; a ^= rot(c,16); a &= 0xffffffff; c += b; c &= 0xffffffff
    b -= a; b &= 0xffffffff; b ^= rot(a,19); b &= 0xffffffff; a += c; a &= 0xffffffff
    c -= b; c &= 0xffffffff; c ^= rot(b,4);  c &= 0xffffffff; b += a; b &= 0xffffffff
    return a, b, c

def final(a, b, c):
    a &= 0xffffffff; b &= 0xffffffff; c &= 0xffffffff
    c ^= b; c &= 0xffffffff; c -= rot(b,14); c &= 0xffffffff
    a ^= c; a &= 0xffffffff; a -= rot(c,11); a &= 0xffffffff
    b ^= a; b &= 0xffffffff; b -= rot(a,25); b &= 0xffffffff
    c ^= b; c &= 0xffffffff; c -= rot(b,16); c &= 0xffffffff
    a ^= c; a &= 0xffffffff; a -= rot(c,4);  a &= 0xffffffff
    b ^= a; b &= 0xffffffff; b -= rot(a,14); b &= 0xffffffff
    c ^= b; c &= 0xffffffff; c -= rot(b,24); c &= 0xffffffff
    return a, b, c

def hashlittle2(data, initval = 0, initval2 = 0):
    def get_byte(a):
        return ord(a)
    if isinstance(data, bytes):
        def get_byte(a):
            return a
    length = lenpos = len(data)

    a = b = c = (0xdeadbeef + (length) + initval)

    c += initval2; c &= 0xffffffff

    p = 0  # string offset
    while lenpos > 12:
        a += (get_byte(data[p+0]) + (get_byte(data[p+1])<<8) + (get_byte(data[p+2])<<16) + (get_byte(data[p+3])<<24)); a &= 0xffffffff
        b += (get_byte(data[p+4]) + (get_byte(data[p+5])<<8) + (get_byte(data[p+6])<<16) + (get_byte(data[p+7])<<24)); b &= 0xffffffff
        c += (get_byte(data[p+8]) + (get_byte(data[p+9])<<8) + (get_byte(data[p+10])<<16) + (get_byte(data[p+11])<<24)); c &= 0xffffffff
        a, b, c = mix(a, b, c)
        p += 12
        lenpos -= 12

    if lenpos == 12: c += (get_byte(data[p+8]) + (get_byte(data[p+9])<<8) + (get_byte(data[p+10])<<16) + (get_byte(data[p+11])<<24)); b += (get_byte(data[p+4]) + (get_byte(data[p+5])<<8) + (get_byte(data[p+6])<<16) + (get_byte(data[p+7])<<24)); a += (get_byte(data[p+0]) + (get_byte(data[p+1])<<8) + (get_byte(data[p+2])<<16) + (get_byte(data[p+3])<<24));
    if lenpos == 11: c += (get_byte(data[p+8]) + (get_byte(data[p+9])<<8) + (get_byte(data[p+10])<<16)); b += (get_byte(data[p+4]) + (get_byte(data[p+5])<<8) + (get_byte(data[p+6])<<16) + (get_byte(data[p+7])<<24)); a += (get_byte(data[p+0]) + (get_byte(data[p+1])<<8) + (get_byte(data[p+2])<<16) + (get_byte(data[p+3])<<24));
    if lenpos == 10: c += (get_byte(data[p+8]) + (get_byte(data[p+9])<<8)); b += (get_byte(data[p+4]) + (get_byte(data[p+5])<<8) + (get_byte(data[p+6])<<16) + (get_byte(data[p+7])<<24)); a += (get_byte(data[p+0]) + (get_byte(data[p+1])<<8) + (get_byte(data[p+2])<<16) + (get_byte(data[p+3])<<24));
    if lenpos == 9:  c += (get_byte(data[p+8])); b += (get_byte(data[p+4]) + (get_byte(data[p+5])<<8) + (get_byte(data[p+6])<<16) + (get_byte(data[p+7])<<24)); a += (get_byte(data[p+0]) + (get_byte(data[p+1])<<8) + (get_byte(data[p+2])<<16) + (get_byte(data[p+3])<<24));
    if lenpos == 8:  b += (get_byte(data[p+4]) + (get_byte(data[p+5])<<8) + (get_byte(data[p+6])<<16) + (get_byte(data[p+7])<<24)); a += (get_byte(data[p+0]) + (get_byte(data[p+1])<<8) + (get_byte(data[p+2])<<16) + (get_byte(data[p+3])<<24));
    if lenpos == 7:  b += (get_byte(data[p+4]) + (get_byte(data[p+5])<<8) + (get_byte(data[p+6])<<16)); a += (get_byte(data[p+0]) + (get_byte(data[p+1])<<8) + (get_byte(data[p+2])<<16) + (get_byte(data[p+3])<<24));
    if lenpos == 6:  b += ((get_byte(data[p+5])<<8) + get_byte(data[p+4])); a += (get_byte(data[p+0]) + (get_byte(data[p+1])<<8) + (get_byte(data[p+2])<<16) + (get_byte(data[p+3])<<24))
    if lenpos == 5:  b += (get_byte(data[p+4])); a += (get_byte(data[p+0]) + (get_byte(data[p+1])<<8) + (get_byte(data[p+2])<<16) + (get_byte(data[p+3])<<24));
    if lenpos == 4:  a += (get_byte(data[p+0]) + (get_byte(data[p+1])<<8) + (get_byte(data[p+2])<<16) + (get_byte(data[p+3])<<24))
    if lenpos == 3:  a += (get_byte(data[p+0]) + (get_byte(data[p+1])<<8) + (get_byte(data[p+2])<<16))
    if lenpos == 2:  a += (get_byte(data[p+0]) + (get_byte(data[p+1])<<8))
    if lenpos == 1:  a += get_byte(data[p+0])
    a &= 0xffffffff; b &= 0xffffffff; c &= 0xffffffff
    if lenpos == 0: return c, b

    a, b, c = final(a, b, c)

    return c, b

def lookup3(data, initval=0):
    c, b = hashlittle2(data, initval, 0)
    return c


if __name__ == "__main__":
    for subdir, dirs, files in os.walk("."):
        for f in files:
            if f.endswith("_master.h5"):
                describe(os.path.join(subdir, f).replace("./", ""))