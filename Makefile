# Copyright (C) 2013 DNAnexus, Inc.
#
# This file is part of gatk_recalibrate_and_realign (DNAnexus platform app).
#
# (The MIT Expat License)
#
#   Permission is hereby granted, free of charge, to any person obtaining a
#   copy of this software and associated documentation files (the "Software"),
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense,
#   and/or sell copies of the Software, and to permit persons to whom the
#   Software is furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in
#   all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.

all: bamtools samtools

bamtools:
	mkdir -p bamtools/build
	cd bamtools/build; cmake ..; make
	cp -a bamtools/bin/bamtools-* resources/usr/bin

samtools:
	make -C samtools samtools
	cp -a samtools/samtools resources/usr/bin/samtools

clean:
	rm -f resources/usr/bin/samtools resources/usr/bin/bamtools*

.PHONY: all bamtools samtools clean
