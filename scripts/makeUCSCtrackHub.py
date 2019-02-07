#! /usr/bin/env python

__author__ = 'AF'

# check the tutorial here https://pythonhosted.org/trackhub/tutorial.html
# import the components we'll be using
# 01/12/2016

import argparse

import glob, os
import trackhub

parser = argparse.ArgumentParser()
parser.add_argument("--hub_name", help="Required. the name of your track hub")
parser.add_argument("--bw", nargs='+', help = "BigWig files to be added to the hub")
parser.add_argument("--peaks", nargs = '+', help = "Called peaks in bed format")
parser.add_argument("--categories", nargs='+', help= "TF or Histone marks. If provided, file names will be searched for pattern provided by this option. Then tracks will be aggregated based on those categories")
parser.add_argument("--output_dir", help="the folder where the hub track files are generated, default is the same as input_dir", default=".")
#parser.add_argument("--email", help="Required. your email to contact")
#parser.add_argument("--composite_track_name", help="Required. the name of your composite track")

args = parser.parse_args()

#assert args.hub_name is not None, "please provide the hub_name"
#assert args.base_url is not None, "please provide the base_url"
#assert args.composite_track_name is not None, "please provide the composite track name"
#assert args.email is not None, "please provide your email"
#assert args.input_dir is not None, "please provide the path to the bigwig and bigbed files on your local computer"

def subgroups_from_filename(fn):
    """
    This functions figures out subgroups based on the number in the
    filename.  Subgroups provided to the Track() constructor is
    a dictionary where keys are `name` attributes from the subgroups added
    to the composite above, and values are keys of the `mapping` attribute
    of that same subgroup.

    Might be easier to cross-reference with the subgroups above, but an
    example return value from this function would be::

        {'s1': 'n', 'num': '2'}
    """
    if fn.endswith('bw'):
        kind = 'signal'
    else:
        kind = 'peak'
    track_subgroup = {
        'File type': kind
    }
    return track_subgroup

# First we initialize the components of a track hub
hub, genomes_file, genome, trackdb = trackhub.default_hub(
    hub_name=args.hub_name,
    short_label=args.hub_name,
    long_label=args.hub_name,
    genome="hg19",
    email="dalerr@niddk.nih.gov")

##Subgroup definition
subgroups = [

    trackhub.SubGroupDefinition(
        name='File type',
        label='File type',
        mapping={
            'signal': 'signal',
            'peak': 'peak',
        }
    )
]

#####Composite track using the different subgroup
composite = trackhub.CompositeTrack(
    name='Composite',
    short_label='Signal and regions',
    # The available options for dimensions are the `name` attributes of
    # each subgroup. Start with dimX and dimY (which become axes of the
    # checkbox matrix to select tracks), and then dimA, dimB, etc.
    dimensions='dimA=kind',

    # This enables a drop-down box under the checkbox matrix that lets us
    # select whatever dimA is (here, "kind").
    filterComposite='dimA',

    # The availalbe options here are the `name` attributes of each subgroup.
    tracktype='bigWig',
    visibility='full',
)
composite.add_subgroups(subgroups)
trackdb.add_tracks(composite)
# CompositeTracks compose different ViewTracks. We'll make one ViewTrack
# for signal, and one for bigBed regions.
signal_view = trackhub.ViewTrack(
    name='signalviewtrack',
    view='signal',
    visibility='full',
    tracktype='bigWig',
    short_label='Signal')

if args.peaks is not None:
    regions_view = trackhub.ViewTrack(
        name='regionsviewtrack',
        view='regions',
        visibility='dense',
        tracktype='bigBed',
        short_label='Regions')
    composite.add_view(regions_view)


# These need to be added to the composite.
composite.add_view(signal_view)

for bigwig in args.bw:
    track = trackhub.Track(
        name=trackhub.helpers.sanitize(os.path.basename(bigwig)),
        source=bigwig,
        visibility='full',
        tracktype='bigWig',
        viewLimits='-2:2',
        maxHeightPixels='8:50:128',
        subgroups=subgroups_from_filename(bigwig),
    )
    # Note that we add the track to the *view* rather than the trackDb as
    # we did in the README example.
    signal_view.add_tracks(track)

# Same thing with the bigBeds. No overlay track to add these to, though.
# Just to the regions_view ViewTrack.
if args.peaks is not None:
    for bigbed in args.peaks:
        track = trackhub.Track(
            name=trackhub.helpers.sanitize(os.path.basename(bigbed)),
            source=bigbed,
            tracktype='bigBed',
            visibility='dense',
            subgroups=subgroups_from_filename(bigbed),
            )
        regions_view.add_tracks(track)

#####Super traxk with Aggregate track
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

                name = trackhub.helpers.sanitize(os.path.basename(bigwig)) + 'agg'
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