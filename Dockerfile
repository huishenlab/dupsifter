################################################################################
##
## THE MIT LICENSE
##
## Copyright 2023 Jacob Morrison <jacob.morrison@vai.org>
##
## Permission is hereby granted, free of charge, to any person obtaining a copy
## of this software and associated documentation files (the “Software”), to deal
## in the Software without restriction, including without limitation the rights
## to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
## copies of the Software, and to permit persons to whom the Software is
## furnished to do so, subject to the following conditions:
##
## The above copyright notice and this permission notice shall be included in
## all copies or substantial portions of the Software.
##
## THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
## IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
## FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
## AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
## LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
## OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
## SOFTWARE.
##
## DESCRIPTION:
##     Container file for building and installing dupsifter
##
## CREATED BY:
##     Jacob Morrison
##
## CREATION DATE:
##     Sep 2023
##
## UPDATE NOTES:
##     - Sep 2023
##         - Initial creation
##
################################################################################

# Base image
FROM ubuntu:latest

# Author/maintainer
MAINTAINER Jacob Morrison

# Set up dependencies
RUN apt update --fix-missing && \
    apt install gcc make curl \
    libcurl4-openssl-dev zlib1g-dev libbz2-dev liblzma-dev

RUN curl -OL \
    $(curl -s https://api.github.com/repos/huishenlab/dupsifter/releases/latest | \
        grep browser_download_url | grep release-source.zip | cut -d '"' -f 4) && \
    unzip release-source.zip && \
    cd dupsifter && \
    make
