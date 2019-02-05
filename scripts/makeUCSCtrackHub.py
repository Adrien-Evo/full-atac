#! /usr/bin/env python

__author__ = 'tommy'

# check the tutorial here https://pythonhosted.org/trackhub/tutorial.html
# import the components we'll be using
# 01/12/2016

import argparse

import glob, os
import trackhub

parser = argparse.ArgumentParser()
parser.add_argument("--hub_name", help="Required. the name of your track hub")
parser.add_argument("--bw", nargs='+', help=" BigWig files to be added to the hub")
parser.add_argument("--categories", nargs='+', help= "If provided, file names will be searched for pattern provided by this option. Then tracks will be aggregated based on those categories")
parser.add_argument("--output_dir", help="the folder where the hub track files are generated, default is the same as input_dir", default=".")
#parser.add_argument("--email", help="Required. your email to contact")
#parser.add_argument("--composite_track_name", help="Required. the name of your composite track")

args = parser.parse_args()

#assert args.hub_name is not None, "please provide the hub_name"
#assert args.base_url is not None, "please provide the base_url"
#assert args.composite_track_name is not None, "please provide the composite track name"
#assert args.email is not None, "please provide your email"
#assert args.input_dir is not None, "please provide the path to the bigwig and bigbed files on your local computer"

# First we initialize the components of a track hub
hub, genomes_file, genome, trackdb = trackhub.default_hub(
    hub_name=args.hub_name,
    short_label=args.hub_name,
    long_label=args.hub_name,
    genome="hg19",
    email="dalerr@niddk.nih.gov")

if args.categories is not None :
#If there is more than 1 category, create a super track to hold all aggregate tracks
    if len(args.categories) > 1 :
        supertrack = trackhub.SuperTrack(
        name = 'SuperTrack',
        short_label = 'SuperTrack'
        )
        trackdb.add_tracks(supertrack)

##Create subgroup from categories. Create a super track with one overlay track per categories of marks or TF
    
    for category in args.categories:
        
        overlay = trackhub.AggregateTrack(
        aggregate='transparentOverlay',
        visibility='full',
        tracktype='bigWig',
        viewLimits='-10:10',
        maxHeightPixels='8:80:128',
        showSubtrackColorOnUi='on',
        name=category)

        if len(args.categories) > 1 :
            supertrack.add_tracks(overlay)
        elif len(args.categories) == 1: 
            trackdb.add_tracks(overlay)
        
        #Now let slooks at the big wig files added and look for the category pattern and add them to the overlay track
        #category string uses the "_" from the main snakemake file that defines samples based on {sample}_{marks}.
        category_string = "_"+category +"_"

        for bigwig in args.bw :
            if category_string in bigwig:

                name = trackhub.helpers.sanitize(os.path.basename(bigwig))
                track = trackhub.Track(
                    name=name,          # track names can't have any spaces or special chars.
                    source= bigwig,      # filename to build this track from
                    visibility='full',  # shows the full signal
                    color='128,0,5',    # brick red
                    autoScale='on',     # allow the track to autoscale
                    tracktype='bigWig', # required when making a track
                )
                
                overlay.add_subtrack(track)
        # Each track is added to the trackdb



print(trackdb)
# In this example we "upload" the hub locally. Files are created in the
# "example_hub" directory, along with symlinks to the tracks' data files.
# This directory can then be pushed to GitHub or rsynced to a server.
#trackhub.upload.upload_hub(hub=hub,host="localhost",remote_dir=args.output_dir)
trackhub.upload.stage_hub(hub=hub,staging=args.output_dir)


# Alternatively, we could upload directly to a web server (not run in this
# example):
#if 0:
#    trackhub.upload.upload_hub(
#        hub=hub, host='example.com', user='username',
#        remote_dir='/var/www/example_hub')