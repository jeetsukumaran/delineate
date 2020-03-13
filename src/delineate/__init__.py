#! /usr/bin/env python

##############################################################################
## Copyright (c) 2018 Jeet Sukumaran and Mark T. Holder
## All rights reserved.
##
## Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are met:
##
##     * Redistributions of source code must retain the above copyright
##       notice, this list of conditions and the following disclaimer.
##     * Redistributions in binary form must reproduce the above copyright
##       notice, this list of conditions and the following disclaimer in the
##       documentation and/or other materials provided with the distribution.
##     * The names of its contributors may not be used to endorse or promote
##       products derived from this software without specific prior written
##       permission.
##
## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
## ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
## WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
## DISCLAIMED. IN NO EVENT SHALL JEET SUKUMARAN BE LIABLE FOR ANY DIRECT,
## INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
## BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
## DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
## LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
## OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
## ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
##
##############################################################################

__project__ = "DELINEATE"
__version__ = "1.2.2"

def revision_description():
    import os
    import subprocess
    import sys
    cwd = os.path.dirname(os.path.abspath(__file__))
    try:
        desc = subprocess.check_output(
                ["git", "log", "-1", "--pretty=format:%h, %ci"],
                stderr=subprocess.PIPE,
                cwd=cwd).strip().decode(sys.stdout.encoding)[:-6]
    except:
        return ""
    if desc:
        revision_text = " (revision: {})".format(desc)
    else:
        revision_text = ""
    return revision_text

def name():
    return "{} {}{}".format(__project__, __version__, revision_description())

def homedir():
    import os
    try:
         __homedir__ = os.path.dirname(os.path.abspath(__file__))
    except:
        __homedir__ = None
    return __homedir__

def description(dest=None):
    import sys
    import site
    if dest is None:
        dest = sys.stdout
    fields = {}
    fields["{} version".format(__project__)] = name()
    fields["{} location".format(__project__)] = homedir()
    fields["Python version"] = sys.version.replace("\n", "")
    fields["Python executable"] = sys.executable
    try:
        fields["Python site packages"] = site.getsitepackages()
    except:
        pass
    max_fieldname_len = max(len(fieldname) for fieldname in fields)
    for fieldname, fieldvalue in fields.items():
        dest.write("{fieldname:{fieldnamewidth}}: {fieldvalue}\n".format(
            fieldname=fieldname,
            fieldnamewidth=max_fieldname_len + 2,
            fieldvalue=fieldvalue))

def get_metadata_parser_opts():
    # To use in main(), instantiate the parser by:
    #
    #       parser = argparse.ArgumentParser(description=__doc__,
    #                   parents=[delineate.get_metadata_parser_opts()])
    import argparse
    metadata_parser_opts = argparse.ArgumentParser(add_help=False)
    metadata_parser_opts.add_argument("--version",
            action="store_true",
            default=False,
            help="Print version and exit.")
    metadata_parser_opts.add_argument("--describe",
            action="store_true",
            default=False,
            help="Print package description and exit.")
    return metadata_parser_opts


