#!/bin/env python

from pycirclize import Circos
from pycirclize.parser import Genbank
import numpy as np
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import glob
from Bio import SeqIO
import os 
import argparse
import matplotlib.pyplot as plt

# creating the circulize class
class circluar_genome():
    def __init__(self, gbk_path, space, dpi):
        self.gbk_path = gbk_path
        self.space = space
        self.dpi = dpi
        #self.vc_path = vc_path
        #creating the initial layout
        self.gbk_ref = Genbank(self.gbk_path)
        self.circos = Circos(sectors={self.gbk_ref.name: self.gbk_ref.genome_length + 10000}, space = self.space)
        self.sector = self.circos.get_sector(self.gbk_ref.name)
        
    def plot_layout(self, text: str, radious = 20, tick_numbers = 10, sector_size = (98,100), fill_color = "lightgrey"):
        
        self.circos.text(text, r = radious)
        #plot outer track with xticks
        tic_interval = self.gbk_ref.genome_length / tick_numbers
        outr_track = self.sector.add_track(sector_size)
        outr_track.axis(fc= fill_color)
        outr_track.xticks_by_interval(
            tic_interval, label_formatter=lambda v: f"{v/10**6:.1f} Mb",
            show_label=True
        )
        
        #self.circos.plotfig(dpi=400)
    # adding the cds 
    def cds_add(self, gbk = None, sector_size = (95,97), r_pad_ratio= 0.1, name="CDS", color = "black", 
                strand= 1, label = "CDS", label_color = "black", label_size = 3):
        if gbk is None:
           gbk_obj = Genbank(self.gbk_path)
        else:
            gbk_obj = Genbank(gbk)
        
        if strand ==1:
            f_cds_track = self.sector.add_track(sector_size, r_pad_ratio=r_pad_ratio)
            f_cds_track.genomic_features(gbk_obj.extract_features(feature_type=name, target_strand=strand), fc = color)
            if label is not None:
                self.circos.text( text = label, r = sum(sector_size)/2, size = label_size, color = label_color,  horizontalalignment =  'right')
        elif strand == -1:
            r_cds_track = self.sector.add_track(sector_size, r_pad_ratio=r_pad_ratio)
            r_cds_track.genomic_features(gbk_obj.extract_features(name, target_strand=strand), fc = color)
            if label is not None:
                self.circos.text( text = label, r = sum(sector_size)/2, size = label_size, color = label_color,  horizontalalignment =  'right')
        else:
            raise ValueError("Strand must be either +1 or -1")
        
        #self.circos.plotfig(dpi=self.dpi)
    
    #Adding GC content 
    def gc_add(self, gbk = None, sector_size = (95,97), r_pad_ratio= 0.1, gc_negative_color = "darkgrey", gc_positive_color = "forestgreen", label = "GC content", label_color = "black", label_size = 3):
        if gbk is None:
           gbk_obj = Genbank(self.gbk_path)
        else:
            gbk_obj = Genbank(gbk)

        gc_content_track = self.sector.add_track(sector_size, r_pad_ratio=r_pad_ratio)
        pos_list, gc_contents = gbk_obj.calc_gc_content()
        gc_contents = gc_contents - gbk_obj.calc_genome_gc_content()
        positive_gc_contents = np.where(gc_contents > 0, gc_contents, 0)
        negative_gc_contents = np.where(gc_contents < 0, gc_contents, 0)
        abs_max_gc_content = np.max(np.abs(gc_contents))
        vmin, vmax = -abs_max_gc_content, abs_max_gc_content
        gc_content_track.fill_between(
            pos_list, positive_gc_contents, 0, vmin=vmin, vmax=vmax, color = gc_positive_color
        )
        gc_content_track.fill_between(
            pos_list, negative_gc_contents, 0, vmin=vmin, vmax=vmax, color = gc_negative_color
        )
        self.circos.text( text = label, r = sum(sector_size)/2, size = label_size, color = label_color, horizontalalignment =  'right')
    
    #Adding GC skew
    def gc_skew_add(self, gbk = None, sector_size = (95,97), r_pad_ratio= 0.1,  positive_color = "olive", negative_color = "darkgrey",
        label = "GC skew", label_color = "black", label_size = 3):

        if gbk is None:
           gbk_obj = Genbank(self.gbk_path)
        else:
            gbk_obj = Genbank(gbk)
        gc_skew_track = self.sector.add_track(sector_size, r_pad_ratio=r_pad_ratio)
        pos_list, gc_skews = gbk_obj.calc_gc_skew()
        positive_gc_skew = np.where(gc_skews > 0, gc_skews, 0)
        negative_gc_skew = np.where(gc_skews < 0, gc_skews, 0)
        abs_max_gc_skew = np.max(np.abs(gc_skews))
        vmin, vmax = -abs_max_gc_skew, abs_max_gc_skew
        gc_skew_track.fill_between(
            pos_list, positive_gc_skew, 0, vmin=vmin, vmax=vmax, color = positive_color
        )
        gc_skew_track.fill_between(
            pos_list, negative_gc_skew, 0, vmin=vmin, vmax=vmax, color = negative_color
        )
        self.circos.text( text = label, r = sum(sector_size)/2, size = label_size, color = label_color,   horizontalalignment =  'right')
        
    #adding the variant calls 
    def vc_add(self, vc_df: str, sector_size: (93,95), r_pad_ratio = 0.1, ref_color = "darkmagenta", vc_lab_size = 1.3, vc_color = "olive", ref_catpion = "REF", alt_caption = "ALT", label_size = 3):
        vc_track = self.sector.add_track(sector_size, r_pad_ratio=r_pad_ratio)
        pos_list, ref = vc_df["POS"], vc_df["REF"]
        vc_track.xticks(x = list(pos_list), labels= list(ref), label_size=vc_lab_size, 
                        tick_length=0.60, label_margin = 0.1, label_orientation="horizontal", 
                        text_kws=dict(color = ref_color))
        
        self.circos.text( text = ref_catpion, r =  sector_size[1]+0.9, size = label_size, color = ref_color, horizontalalignment =  'right')
        self.circos.text( text = alt_caption, r = sector_size[0]+0.5, size = label_size, color = vc_color,  horizontalalignment =  'right')

        vc_track = self.sector.add_track(tuple([v + 0.7 for v in sector_size]), r_pad_ratio=r_pad_ratio)
        pos_list, alt = vc_df["POS"], vc_df["ALT"]
        vc_track.xticks(x = list(pos_list), outer=False, labels= list(alt), 
                        label_size=vc_lab_size, tick_length=0, label_orientation="horizontal", 
                        label_margin=0, text_kws=dict(color = vc_color))
        #self.circos.plotfig(dpi=self.dpi)
    
    def save_plot(self, path, dpi = 500):
        self.circos.savefig(path, dpi=dpi)

def plotter():
    parser = argparse.ArgumentParser(description="Generate a circular genome plot with GC content and skew.")
    parser.add_argument('-d', '--gbk_dir', type= str, required=True, help= "Directory containing annotation folders from Prokka or Bakta which each contains a gbk file.")
    parser.add_argument('-o', '--out_dir', type= str, required=True, help= "Output directory.")
    parser.add_argument("--add_gc", default=True, action="store_true", help="Include GC content in the plot.")
    parser.add_argument("--add_skew", action="store_true", default=True, help="Include GC skew in the plot.")
    parser.add_argument("--dpi", type=int, default=300, help="Size of the figure.")
    parser.add_argument("--interval", type=int, default=3, help="Interval for sector size adjustment.")
    parser.add_argument("--figsize", type = int, default=10, help="Figure size (inc)")
    parser.add_argument("--f_color", type=str, default="deepskyblue", help="Forward strand color")
    parser.add_argument("--r_color", type=str, default="#ff7261", help="Reverse strand color")

    args = parser.parse_args()

    files = glob.glob(os.path.join(args.gbk_dir, "*/*"))
    col_forward = args.f_color
    col_reverse = args.r_color
    out_file = os.path.join(args.out_dir, "circular_plot.png")
    add_gc = args.add_gc
    add_skew = args.add_skew
    dpi = args.dpi
    interval = args.interval

    gbks = []
    lenghts = {}
    for file in files:
        if file.endswith(".gbk"):
            gbks.append(file)
            recs = list(SeqIO.parse(file, 'genbank'))
            leng = sum(len(record.seq) for record in recs)# get total length sum
            lenghts[file] = leng
    start_genome = max(lenghts, key=lenghts.get)

    pl = circluar_genome(start_genome, dpi=dpi, space=40)
    pl.plot_layout(
        text=None 
        )


    ref_name  = os.path.splitext(os.path.basename(start_genome))[0]
    pl.cds_add(
        strand = 1,  
        sector_size= ((97-interval),97), 
        color= col_forward, 
        label=ref_name, 
        label_size = 7
        )
    pl.cds_add(
        strand = -1, 
        sector_size= ((97-interval),97),  
        color= col_reverse, 
        label=None)


    num_file = len(gbks) - 1

    ranges = []
    start = 97 - interval 

    for _ in range(num_file):
        end = start
        start -=  interval
        ranges.append((start, end))



    #remove the reference
    rm = gbks.index(start_genome)

    gbks.pop(rm)


    for file, index in zip(gbks, range(0, len(gbks))):
        sector_name = os.path.splitext(os.path.basename(file))[0]
        name = sector_name
        
        pl.cds_add(gbk=  file, strand = 1, sector_size=ranges[index], color= col_forward, label = name, label_size= 7)
        pl.cds_add(gbk=  file, strand = -1, sector_size=ranges[index], color=col_reverse, label=None)



    #Adding GC

    if add_gc:
        gc_tuple = []

        start = ranges[-1][-1] - interval

        for _ in range(num_file):
            end = start 
            start -= interval
            gc_tuple.append((start, end)) 

        for file, index in zip(gbks, range(0, len(gbks))):
            sector_name = os.path.splitext(os.path.basename(file))[0]
            name = sector_name

            pl.gc_add(label=f"GC content: {name}", sector_size=gc_tuple[index], label_size= 5)

        

    # gc skew
    if add_skew:
        skew_tuple = []

        start = gc_tuple[-1][-1] - interval

        for _ in range(num_file):
            end = start 
            start -= interval
            skew_tuple.append((start, end)) 

        for file, index in zip(gbks, range(0, len(gbks))):
            sector_name = os.path.splitext(os.path.basename(file))[0]
            name = sector_name

            pl.gc_skew_add(label=f"GC skew: {name}", sector_size=skew_tuple[index], label_size= 5)



    pl.circos.plotfig(dpi=dpi, figsize=(args.figsize, args.figsize ))

    if add_gc and add_skew:
        handles = [
            Patch(color = "white", label =  f"Pangenome of {len(lenghts)} genomes"),
            Patch(color = col_forward, label = "Forward CDS"),
            Patch(color=col_reverse, label = "Reverse CDS"),
            Line2D([], [], color="forestgreen", label = "Positive GC content", marker="^", ms = 5, ls = "None"),
            Line2D([], [], color="grey", label="Negative GC Content", marker="v", ms=5, ls="None"),
            Line2D([], [], color="olive", label="Positive GC Skew", marker="^", ms=5, ls="None"),
            Line2D([], [], color="grey", label="Negative GC Skew", marker="v", ms=5, ls="None")

        ]
    elif add_gc and not add_skew:
        handles = [
            Patch(color = "white", label =  f"Pangenome of {len(lenghts)} genomes"),
            Patch(color = col_forward, label = "Forward CDS"),
            Patch(color=col_reverse, label = "Reverse CDS"),
            Line2D([], [], color="forestgreen", label = "Positive GC content", marker="^", ms = 5, ls = "None"),
            Line2D([], [], color="grey", label="Negative GC Content", marker="v", ms=5, ls="None")

        ]
    elif add_skew and not add_gc:
        handles = [
            Patch(color = "white", label =  f"Pangenome of {len(lenghts)} genomes"),
            Patch(color = col_forward, label = "Forward CDS"),
            Patch(color=col_reverse, label = "Reverse CDS"),
            Line2D([], [], color="olive", label="Positive GC Skew", marker="^", ms=5, ls="None"),
            Line2D([], [], color="grey", label="Negative GC Skew", marker="v", ms=5, ls="None")

        ]
    else:
        handles = [
            Patch(color = "white", label =  f"Pangenome of {len(lenghts)} genomes"),
            Patch(color = col_forward, label = "Forward CDS"),
            Patch(color=col_reverse, label = "Reverse CDS")

        ]


    pl.circos.ax.legend(
    handles=handles, 
        bbox_to_anchor = (0.5, 0.5), 
        loc = "center", 
        fontsize = 10
        ).figure.savefig(out_file, dpi = 500);
    print(f"✅ Circular plot saved at {args.out_dir}{out_file}")
if __name__ == "__main__":
    plotter()