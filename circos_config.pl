use strict;
use warnings;
use Getopt::Std;
our ($opt_p, $opt_v, $opt_c, $opt_s, $opt_o);
print "RAW paramters: @ARGV\n";
getopt('pvcso');
if ( (!defined $opt_p) && (!defined $opt_v) && (!defined $opt_c) && (!defined $opt_s) && (!defined $opt_o) ) {
    die ("Usage: $0\n\t -p [script_path] \n\t -v[structural vaiant file ] \n\t -c [cnv file] \n\t -s [sample] \n\t -o [outputdir]  \n");
}
else    {
    my $path=$opt_p;
    my $sv_file=$opt_v;
    my $cnv_file=$opt_c;
    my $sample=$opt_s;
    my $output=$opt_o;
    open OUT, ">$output/$sample.main.config" or die "";
    print OUT "<colors>
	<<include ${path}/colors.conf>>
    </colors>
    
    <fonts>
	<<include ${path}/fonts.conf>>
    </fonts>
    
    <<include ${path}/fusion_ideogram.conf>>
    <<include ${path}/ticks.conf>>
    
    karyotype = ${path}/karyotype.conf
    
    <links>
    
    z = 0
    radius = 0.73r
    bezier_radius = 0.2r
    
    <link mylink>
    show = yes
    color = grey
    thickness = 1
    file = ${sv_file}
    
    <rules>
    <rule>
    importance = 90
    condition = _INTRACHR_
    color = red
    z = 90
    </rule>
    <rule>
    condition = _INTERCHR_
    color = blue
    bezier_radius = 0.05r
    </rule>
    </rules>
    </link>
    
    </links>
    
    
    <plots>
    <plot>
    
    show    = yes
    type    = line
    max_gap = 10u
    file    = ${cnv_file}
    color = vdgrey
    thickness = 1
    min   = 0
    max   = 2.5
    r0    = 0.8r
    r1    = 0.975r
    
    background       = yes
    background_color = vvlgrey
    background_stroke_color = black
    background_stroke_thickness = 1
    
    axis           = yes
    axis_color     = lgrey
    axis_thickness = 2
    axis_spacing   = 0.001
    
    <rules>
    
    <rule>
    importance   = 100
    condition    = _VALUE_ > 1.5
    color = orange
    </rule>
    
    <rule>
    importance   = 85
    condition    = _VALUE_ < 0.5
    color = blue
    </rule>
    
    </rules>
    
    </plot>
    
    </plots>
    
    
    <image>	
	dir = $output
	file = $sample.sv_cnv.png
	png = yes
	24bit = yes
	radius = 1200p
	background = white
	angle_offset = -90
	auto_alpha_colors = yes
	auto_alpha_steps = 5
    </image>
    
    chromosomes_units = 1000000
    chromosomes_display_default = yes
    
    anglestep       = 0.5
    minslicestep    = 10
    beziersamples   = 40
    debug           = no
    warnings        = no
    imagemap        = no
    
    units_ok = bupr
    units_nounit = n";
    close OUT;
}


## end of script
    
