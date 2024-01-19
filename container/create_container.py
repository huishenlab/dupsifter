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
##     Create a Dockerfile for building and installing dupsifter
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
import argparse
import os

LATEST = '1.2.1'

def valid_versions():
    """Valid versions of dupsifter and their corresponding tags.

    Inputs -
        None
    Returns -
        dict
    """
    return {
        '1.0.0': 'v1.0.0.20220804',
        '1.1.0': 'v1.1.0.20230608',
        '1.1.1': 'v1.1.1.20230615',
        '1.2.0': 'v1.2.0.20230926',
        '1.2.1': 'v1.2.1.20240119',
    }

def cli():
    """Command line interface.

    Inputs -
        None
    Returns
        argparse.ArgumentParser
    """
    parser = argparse.ArgumentParser(description='Create Dockerfile for specified input version of dupsifter')

    valid = valid_versions()
    parser.add_argument(
        '-v', '--version',
        type=str,
        help='Version to create Dockerfile for',
        choices=list(valid.keys()), default=LATEST
    )
    parser.add_argument(
        '-f', '--force',
        help='Overwrite existing Dockerfile',
        action='store_true'
    )

    return parser.parse_args()

def create_text(tag):
    """Create Dockerfile using input version.

    Inputs -
        version - tag of dupsifter version
    Returns -
        str
    """
    return '\n'.join(
        [
            '# Base image',
            'FROM ubuntu:latest',
            '',
            '# Set up dependencies',
            'RUN apt update --fix-missing && \\',
            '    apt install -qy gcc make wget unzip \\',
            '    libcurl4-openssl-dev zlib1g-dev libbz2-dev liblzma-dev',
            '',
            '# Build dupsifter',
            f'RUN wget https://github.com/huishenlab/dupsifter/releases/download/{tag}/release-source.zip && \\',
            '    unzip release-source.zip && \\',
            '    cd dupsifter-release && \\',
            '    make',
            '',
            '# Add dupsifter to PATH',
            'ENV PATH dupsifter-release:$PATH',
        ]
    )

def main():
    """Create Dockerfile.

    Inputs -
        None
    Returns -
        None
    """
    args = cli()

    fname = f'Dockerfile_v{args.version}'
    print(f'Will attempt to create {fname} ...')
    if not args.force and os.path.exists(fname):
        print(f'ERROR: {fname} exists already. Please either rename or use -f/--force to overwrite.')
        return None

    valid = valid_versions()
    tag = valid[args.version]

    output = create_text(tag)
    with open(fname, 'w') as fh:
        fh.write(f'{output}')
    print(f'{fname} created successfully')

    return None

if __name__ == '__main__':
    main()
