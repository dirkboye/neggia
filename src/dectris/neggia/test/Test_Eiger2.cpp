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

#include <dectris/neggia/user/Dataset.h>
#include <dectris/neggia/user/H5File.h>
#include "DatasetsFixture.h"

TEST_F(TestDatasetEiger2001, KeepsFileOpen) {
    Dataset xp(H5File(getPathToSourceFile()),
               "/entry/instrument/detector/x_pixel_size");
    double val;
    xp.read(&val);
    ASSERT_EQ(val, double(7.5e-05));
}

TEST_F(TestDatasetEiger2001, MasterFile) {
    {
        Dataset xp(H5File(getPathToSourceFile()),
                   "/entry/instrument/detector/x_pixel_size");
        assert(xp.dim().empty());
        assert(xp.isChunked() == false);
        assert(xp.dataTypeId() == 1);
        assert(xp.dataSize() == sizeof(double));
        double val;
        xp.read(&val);
        ASSERT_EQ(val, double(7.5e-05));
    }
    {
        Dataset yp(H5File(getPathToSourceFile()),
                   "/entry/instrument/detector/y_pixel_size");
        assert(yp.dim().empty());
        assert(yp.isChunked() == false);
        assert(yp.dataTypeId() == 1);
        assert(yp.dataSize() == sizeof(double));
        double val;
        yp.read(&val);
        ASSERT_EQ(val, double(7.5e-05));
    }
    {
        Dataset pixelMask(
                H5File(getPathToSourceFile()),
                "/entry/instrument/detector/detectorSpecific/pixel_mask");
    }
}

TEST_F(TestDatasetEiger2001, DataFile) {
    for (size_t datasetid = 0; datasetid < getNumberOfDatasets(); ++datasetid) {
        Dataset dataset(H5File(getPathToSourceFile()),
                        getTargetDataset(datasetid));
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    ::testing::GTEST_FLAG(catch_exceptions) = false;
    return RUN_ALL_TESTS();
}
