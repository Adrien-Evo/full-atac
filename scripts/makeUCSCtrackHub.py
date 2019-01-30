#! /usr/bin/env python

__author__ = 'tommy'

# check the tutorial here https://pythonhosted.org/trackhub/tutorial.html
# import the components we'll be using
# 01/12/2016

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--hub_name", help="Required. the name of your track hub")
#parser.add_argument("--base_url", help="Required. the path to the folder where the files are accessible from a web, "
#                                       "make sure add the trailing slash")
parser.add_argument("--input_dir", help="Required. the folder where the files are stored in your local computer")
parser.add_argument("--output_dir", help="the folder where the hub track files are generated, default is "
                                         "the same as input_dir", default=".")
#parser.add_argument("--email", help="Required. your email to contact")
#parser.add_argument("--composite_track_name", help="Required. the name of your composite track")

args = parser.parse_args()

assert args.hub_name is not None, "please provide the hub_name"
#assert args.base_url is not None, "please provide the base_url"
#assert args.composite_track_name is not None, "please provide the composite track name"
#assert args.email is not None, "please provide your email"
assert args.input_dir is not None, "please provide the path to the bigwig and bigbed files on your local computer"


import glob, os
import trackhub

# First we initialize the components of a track hub
hub, genomes_file, genome, trackdb = trackhub.default_hub(
    hub_name="test",
    short_label='myhub',
    long_label='myhub',
    genome="hg19",
    email="dalerr@niddk.nih.gov")

# Next, we add a track for every bigwig found.  In practice, you would
# point to your own files. In this example we use the path to the data
# included with trackhub.

for bigwig in glob.glob(os.path.join(args.input_dir,"*bw")):
    name = trackhub.helpers.sanitize(os.path.basename(bigwig))
    track = trackhub.Track(
        name=name,          # track names can't have any spaces or special chars.
        source= os.path.join("https://bird2cluster.univ-nantes.fr/ucsc/Adrien/",os.path.basename(args.output_dir),bigwig),      # filename to build this track from
        visibility='full',  # shows the full signal
        color='128,0,5',    # brick red
        autoScale='on',     # allow the track to autoscale
        tracktype='bigWig', # required when making a track
    )

    # Each track is added to the trackdb
    trackdb.add_tracks(track)

# In this example we "upload" the hub locally. Files are created in the
# "example_hub" directory, along with symlinks to the tracks' data files.
# This directory can then be pushed to GitHub or rsynced to a server.
trackhub.upload.upload_hub(hub=hub,host="localhost",remote_dir=args.output_dir)

# Alternatively, we could upload directly to a web server (not run in this
# example):
#if 0:
#    trackhub.upload.upload_hub(
#        hub=hub, host='example.com', user='username',
#        remote_dir='/var/www/example_hub')