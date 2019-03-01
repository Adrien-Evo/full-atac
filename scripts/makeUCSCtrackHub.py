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
parser.add_argument("--sample_name",nargs='+', help="Required. name of your samples")
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





# First we initialize the components of a track hub
hub, genomes_file, genome, trackdb = trackhub.default_hub(
    hub_name=args.hub_name,
    short_label=args.hub_name,
    long_label=args.hub_name,
    genome="hg19",
    email="adrien.foucal@univ-nantes.fr")

####Definition of subgroup from filename
def subgroups_from_filename(fn):
    """
    This functions figures out subgroups based on the filename and the provided 
    sample_name, marks or TF name. Subgroups provided to the Track() constructor
    is a dictionary where keys are `name` attributes from the subgroups added
    to the composite track, and values are keys of the `mapping` attribute of 
    that same subgroup.

    An example return value from this function would be:
        {'track': 'signal', 'samp': 'DC1077', 'cat' = 'H3K27'}
    """
    track_subgroup = {}
    #track based on extension
    if fn.endswith('bw'.lower())  or fn.endswith('bigWig'.lower()) :
        track_subgroup['track'] = 'signal'
    else:
        track_subgroup['track']= 'peak'

    ##Checking the file name for the provided samp
    for sample_name in args.sample_name:
        if sample_name in fn:
            track_subgroup['samp'] = sample_name

    if args.categories is not None :
        for category in args.categories:
            if category in fn:
                track_subgroup['cat'] = category
    
    return track_subgroup

def simplify_filename(fn):
    """
    This functions simplify the file name to sampleName_Mark and ensurethere's no  trailing characters, based on sample name and mark/TF argumentsfilename and the provided
    """

    simplifiedName=""
    ##Checking the file name for the provided samp
    for sample_name in args.sample_name:
        if sample_name in fn:
            simplifiedName = simplifiedName + sample_name
    
    if args.categories is not None :
        simplifiedName = simplifiedName + "_"
        for category in args.categories:
            #This is not ideal in a context where it's not used in a the pipeline
            category_string = "_"+category +"_"
            if category in fn:
                simplifiedName = simplifiedName +  category
    
    return simplifiedName

# Sample subgroup
dict_sample_name = {args.sample_name[i]: args.sample_name[i] for i in range(len(args.sample_name))}
if args.categories is not None :
    dict_categories = {args.categories[i]: args.categories[i] for i in range(len(args.categories))}


##Subgroup definition
subgroups = [
    #track subgroup
    trackhub.SubGroupDefinition(
        name='track',
        label='File_Type',
        mapping={
            'signal': 'signal',
            'peak': 'peak',
        }
    ),
    trackhub.SubGroupDefinition(
        name='samp',
        label='Sample_Name',
        mapping=dict_sample_name
    )
]

if args.categories is not None :
    subgroups.append(
        trackhub.SubGroupDefinition(
            name='cat',
            label='Categories',
            mapping=dict_categories
        )
    )

#####Composite track using the different subgroup
composite = trackhub.CompositeTrack(
    name='Composite',
    short_label='Signal and regions',
    # The available options for dimensions are the `name` attributes of
    # each subgroup. Start with dimX and dimY (which become axes of the
    # checkbox matrix to select tracks), and then dimA, dimB, etc.
    dimensions='dimX=samp dimY=cat dimA=track',
    # This enables a drop-down box under the checkbox matrix that lets us
    # select whatever dimA is (here, "kind").
    filterComposite='dimA',

    sortOrder='samp=+ cat=- track=+',
    # The available options here are the `name` attributes of each subgroup.
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
        view='peak',
        visibility='dense',
        tracktype='bigWig',
        short_label='Peaks')
    composite.add_view(regions_view)


# These need to be added to the composite.
composite.add_view(signal_view)

##Gogin through all big wigs files and adding them to the composite track
for bigwig in args.bw:
    track = trackhub.Track(
        name=simplify_filename(os.path.basename(bigwig))+"_signal",
        source=bigwig,
        visibility='full',
        tracktype='bigWig',
        autoScale='on',
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
            name=simplify_filename(os.path.basename(bigbed)) + "_peaks",
            source=bigbed,
            tracktype='bigBed',
            visibility='dense',
            autoScale='on',
            subgroups=subgroups_from_filename(bigbed),
            )
        regions_view.add_tracks(track)


#####Super tracks with Aggregate track
if args.categories is not None :
#If there is more than 1 category, create a super track to hold all aggregate tracks
    if len(args.categories) > 1 :
        supertrack = trackhub.SuperTrack(
        name = 'AggregatedCategories',
        short_label = 'Aggregated_Marks_or_TF'
        )
        trackdb.add_tracks(supertrack)

##Create subgroup from categories. Create a super track with one overlay track per categories of marks or TF
    
    for category in args.categories:
        
        overlay = trackhub.AggregateTrack(
        aggregate='transparentOverlay',
        visibility='full',
        tracktype='bigWig',
        maxHeightPixels='8:80:128',
        showSubtrackColorOnUi='on',
        name=category)

        if len(args.categories) > 1 :
            supertrack.add_tracks(overlay)
        elif len(args.categories) == 1: 
            trackdb.add_tracks(overlay)
        
        #Now let slooks at the big wig files added and look for the category pattern and add them to the overlay track
        #category string uses the "_" from the main snakemake file that defines samples based on {sample}_{marks}.
        #This is not ideal in a context where it's not used in a the pipeline
        category_string = "_"+category +"_"
        for bigwig in args.bw :
            if category_string in bigwig:

                name = simplify_filename(os.path.basename(bigwig)) + '_agg'
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



#print(trackdb)
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