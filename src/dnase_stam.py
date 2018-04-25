#!/usr/bin/env python3

from __future__ import print_function

import os
import sys
import re
import json
import random
from collections import OrderedDict, defaultdict

sys.path.append(os.path.join(os.path.dirname(__file__), '../../metadata/utils/'))
from utils import Utils, eprint, AddPath, printt, printWroteNumLines
from exp import Exp

AssayColors = {"DNase": ["6,218,147", "#06DA93"],
               "RNA-seq": ["0,170,0", "", "#00aa00"],
               "RAMPAGE": ["214,66,202", "#D642CA"],
               "H3K4me1": ["255,223,0", "#FFDF00"],
               "H3K4me2": ["255,255,128", "#FFFF80"],
               "H3K4me3": ["255,0,0", "#FF0000"],
               "H3K9ac": ["255,121,3", "#FF7903"],
               "H3K27ac": ["255,205,0", "#FFCD00"],
               "H3K27me3": ["174,175,174", "#AEAFAE"],
               "H3K36me3": ["0,128,0", "#008000"],
               "H3K9me3": ["180,221,228", "#B4DDE4"],
               "Conservation": ["153,153,153", "#999999"],
               "TF ChIP-seq": ["18,98,235", "#1262EB"],
               "CTCF": ["0,176,240", "#00B0F0"]}

ExclusionLabels = [ "all good", "SPOT only", "Zhiping only", "Zhiping and SPOT", "tSNE only", "SPOT and tSNE", "tSNE and Zhiping", "SPOT, tSNE, and Zhiping" ]

SubGroupKeys = ["age", "donor", "label", "assay", "view"]

# from https://www.w3schools.com/colors/colors_shades.asp
COLORS = ["0000CC", "0000FF", "003300", "003333", "003366", "003399", "0033CC", "0033FF", "006600", "006633", "006666", "006699", "0066CC", "0066FF", "009900", "009933", "009966", "009999", "0099CC", "0099FF", "00CC00", "00CC33", "00CC66", "00CC99", "00CCCC", "00CCFF", "00FF00", "00FF33", "00FF66", "00FF99", "00FFCC", "00FFFF", "330000", "330033", "330066", "330099", "3300CC", "3300FF", "333300", "333333", "333366", "333399", "3333CC", "3333FF", "336600", "336633", "336666", "336699", "3366CC", "3366FF", "339900", "339933", "339966", "339999", "3399CC", "3399FF", "33CC00", "33CC33", "33CC66", "33CC99", "33CCCC", "33CCFF", "33FF00", "33FF33", "33FF66", "33FF99", "33FFCC", "33FFFF", "660000", "660033", "660066", "660099", "6600CC", "6600FF", "663300", "663333", "663366", "663399", "6633CC", "6633FF", "666600", "666633", "666666", "666699", "6666CC", "6666FF", "669900", "669933", "669966", "669999", "6699CC", "6699FF", "66CC00", "66CC33", "66CC66", "66CC99", "66CCCC", "66CCFF", "66FF00", "66FF33", "66FF66", "66FF99", "66FFCC", "66FFFF", "990000", "990033", "990066", "990099", "9900CC", "9900FF", "993300", "993333", "993366", "993399", "9933CC", "9933FF", "996600", "996633", "996666", "996699", "9966CC", "9966FF", "999900", "999933", "999966", "999999", "9999CC", "9999FF", "99CC00", "99CC33", "99CC66", "99CC99", "99CCCC", "99CCFF", "99FF00", "99FF33", "99FF66", "99FF99", "99FFCC", "99FFFF", "CC0000", "CC0033", "CC0066", "CC0099", "CC00CC", "CC00FF", "CC3300", "CC3333", "CC3366", "CC3399", "CC33CC", "CC33FF", "CC6600", "CC6633", "CC6666", "CC6699", "CC66CC", "CC66FF", "CC9900", "CC9933", "CC9966", "CC9999", "CC99CC", "CC99FF", "CCCC00", "CCCC33", "CCCC66", "CCCC99", "CCCCCC", "CCCCFF", "CCFF00", "CCFF33", "CCFF66", "CCFF99", "CCFFCC", "CCFFFF", "FF0000", "FF0033", "FF0066", "FF0099", "FF00CC", "FF00FF", "FF3300", "FF3333", "FF3366", "FF3399", "FF33CC", "FF33FF", "FF6600", "FF6633", "FF6666", "FF6699", "FF66CC", "FF66FF", "FF9900", "FF9933", "FF9966", "FF9999", "FF99CC", "FF99FF", "FFCC00", "FFCC33", "FFCC66", "FFCC99", "FFCCCC", "FFCCFF", "FFFF00", "FFFF33", "FFFF66", "FFFF99", "FFFFCC", "FFFFFF"]
COLORS = [tuple(int(h[i:i+2], 16) for i in (0, 2 ,4)) for h in COLORS]
COLORS= [','.join([str(x) for x in c]) for c in COLORS]
random.seed(18124312)
random.shuffle(COLORS)

def viz(state, active):
    if active:
        return state
    return "hide"

def sanitize(s, replChar='_'):
    return re.sub('[^0-9a-zA-Z]+', replChar, s)

def unrollEquals(sUnsorted):
    r = ''
    s = OrderedDict(sorted(sUnsorted.items()))
    for k, v in s.items():
        r += k + '=' + v + ' '
    return r

def getOrUnknown(s):
    if not s:
        return "unknown"
    return s

def makeTrackName(n):
    n = n.replace(" ", "_").replace('(','').replace(')','')
    n = n[:100]
    return n

def makeLongLabel(n):
    return n[:80]

def makeShortLabel(*n):
    return ' '.join([x for x in n if x])[:17]

html_escape_table = {
    "&": "&amp;",
    '"': "&quot;",
    "'": "&apos;",
    ">": "&gt;",
    "<": "&lt;",
    ".": "&#46;",
    ' ': "&#32;"
}

def html_escape(text):
    """Produce entities within text."""
    return "".join(html_escape_table.get(c,c) for c in text)

def colorize(exp):
    if exp.isDNaseSeq():
        return AssayColors["DNase"][0]
    c = "227,184,136"
    if exp.tf in AssayColors:
        c = AssayColors[exp.tf][0]
    if not c:
        if exp.isChipSeqTF():
            if "CTCF" == exp.tf:
                c = AssayColors["CTCF"][0]
            else:
                c = AssayColors["TF ChIP-seq"][0]
    return c


class DetermineTissue:
    # translate tissue name to tissue name
    lookupTissue = {}
    lookupTissue["hg19"] = {}
    lookupTissue["mm10"] = {"small intestine": "intestine",
                            "large intestine": "intestine",
                            "bone element": "bone"}

    # translate biosample term name
    lookupBTN = {}
    fnp = os.path.join(os.path.dirname(__file__), "cellTypeToTissue.hg19.json")
    lookupBTN["hg19"] = json.load(open(fnp))
    fnp = os.path.join(os.path.dirname(__file__), "cellTypeToTissue.mm10.json")
    lookupBTN["mm10"] = json.load(open(fnp))

    @staticmethod
    def TranslateTissue(assembly, exp):
        t = exp.jsondata.get("organ_slims", "")
        if t:
            t = sorted(t)[0]
        else:
            t = ""
        lookup = DetermineTissue.lookupTissue[assembly]
        if t in lookup:
            return lookup[t]
        ct = exp.biosample_term_name
        lookup = DetermineTissue.lookupBTN[assembly]
        if ct in lookup:
            return lookup[ct]
        ct = exp.jsondata.get("biosample_summary", "")
        if ct in lookup:
            return lookup[ct]
        if ct and ct.endswith("erythroid progenitor cells"):
            return "blood"
        if "ENCSR626RVD" == exp.encodeID:
            return "brain"
        if "ENCSR820WLP" == exp.encodeID:
            return "stem cells"
        eprint(assembly, "missing tissiue assignemnt for", exp.encodeID, exp.biosample_term_name)
        return ""

def outputLines(d, indentLevel, extras = {}):
    prefix = '\t' * indentLevel
    for k, v in d.items():
        if v:
            yield prefix + k + " " + str(v) + '\n'
    for k, v in extras.items():
        if v:
            yield prefix + k + " " + str(v) + '\n'
    yield '\n'

class DNaseStam(object):

    def load_newgroup(self):
        self.newgroup = []
        with open("newgroup.tsv", 'r') as f:
            for line in f:
                line = line.strip().split('\t')
                self.newgroup.append(line[0])
    
    def load_exclusion_labels(self):
        self.labels = {}
        with open("exclusion_labels.tsv", 'r') as f:
            for line in f:
                line = line.strip().split('\t')
                if len(line) < 2: continue
                self.labels[line[0]] = ExclusionLabels[int(line[1])]
                
    def loadExps(self):
        printt("loading exps...")
        with open(os.path.join(os.path.dirname(__file__), "hg38-Hotspot-List.txt")) as f:
            lines = [line.rstrip().split() for line in f]
        
        self.exps = []
        for idx, line in enumerate(lines):
            expID = line[0]
            exp = Exp.fromJsonFile(expID)

            t = DetermineTissue.TranslateTissue("hg19", exp).strip()
            if not t:
                raise Exception("missing " + expID)
            
            self.exps.append((expID, exp, t, line[2]))
        printt("loaded", len(self.exps))
            
    def data(self, tissue, exp, fileID, parent, color, tfsx = ""):
        f = None
        for f in exp.files:
            if fileID == f.fileID:
                break
        if fileID != f.fileID:
            print(exp)
            for f in exp.files:
                print(f.fileID)
            raise Exception("missing", fileID)

        p = OrderedDict()
        p["exclusionlabel"] = self.labels[fileID] if fileID in self.labels else ExclusionLabels[0]
        p["track"] = fileID + tfsx
        p["parent"] = parent + (tfsx if "compos" in parent else "")
        p["bigDataUrl"] = f.url + "?proxy=true"
        p["visibility"] = viz("full", True)
        p["type"] = "bigWig"
        p["color"] = color
        p["height"] = "maxHeightPixels 32:12:8"
        p["shortLabel"] = makeShortLabel(exp.biosample_term_name)
        p["longLabel"] = makeLongLabel(fileID + " "+ exp.biosample_summary)
        p["itemRgb"] = "On"
        p["darkerLabels"] = "on"
        return p

    def _doByTissue(self, f):       
        def out(tissue, exp, fileID, tissueColor):
            parent = "compos_" + sanitize(tissue)
            ret = ""
            for line in outputLines(self.data(tissue, exp, fileID, parent, tissueColor), 1):
                ret += line
            return ret

        tissueColors = {}
        byTissue = defaultdict(list)
        for idx, tup in enumerate(sorted(self.exps, key = lambda x: x[0])):
            expID, exp, t, fileID = tup
            if t not in tissueColors:
                tissueColors[t] = COLORS.pop()
            tissueColor = tissueColors[t]
            print(idx, "of", len(self.exps), t, '\t\t',
                  expID, exp.biosample_term_name, exp.age_display)
            byTissue[t].append(out(t, exp, fileID, tissueColor))

        tissueTracks = []
        for t in sorted(byTissue.keys()):
            stanzas = byTissue[t]
            tissueTracks.append("""
track {tn}
parent super_byTissue
compositeTrack on
shortLabel {shortL}
longLabel {longL}
type bigWig 9 +
visibility full
maxHeightPixels 32:12:8
autoScale on
dragAndDrop subTracks
hoverMetadata on
darkerLabels on
""".format(tn = "compos_" + sanitize(t),
           shortL=makeShortLabel(t),
           longL=makeLongLabel(t + " (" + str(len(stanzas)) + " exps)")))

        for tissueTrack in tissueTracks:
            f.write(tissueTrack)
            f.write('\n')
        for t, stanzas in byTissue.items():
            for stanza in sorted(stanzas):
                f.write(stanza)
                f.write('\n')

    def _doByExclusionLabel(self, f):       
        def out(elabel, tissue, exp, fileID, tissueColor):
            parent = "compos_" + sanitize(elabel)
            ret = ""
            for line in outputLines(self.data(tissue, exp, fileID, parent, tissueColor, tfsx = "_el"), 1):
                ret += line
            return ret

        labelColors = {}; toappendlater = {}
        byELabel = defaultdict(list)
        for idx, tup in enumerate(sorted(self.exps, key = lambda x: x[0])):
            expID, exp, t, fileID = tup
            if t not in labelColors:
                labelColors[t] = COLORS.pop()
            labelColor = labelColors[t]
            print(idx, "of", len(self.exps), t, '\t\t',
                  expID, exp.biosample_term_name, exp.age_display)
            if fileID not in self.labels and fileID in self.newgroup:
                toappendlater[fileID] = tup
                continue
            l = self.labels[fileID] if fileID in self.labels else ExclusionLabels[0]
            byELabel[l].append(out(l, t, exp, fileID, labelColor))

        for fileID in self.newgroup:
            if fileID not in toappendlater: continue
            expID, exp, t, fileID = toappendlater[fileID]
            if t not in labelColors:
                labelColors[t] = COLORS.pop()
            labelColor = labelColors[t]
            print(idx, "of", len(self.exps), t, '\t\t',
                  expID, exp.biosample_term_name, exp.age_display)
            l = "new group"
            byELabel[l].append(out(l, t, exp, fileID, labelColor))

        print(str(len(toappendlater)) + " files in new group")
            
        labelTracks = []
        for l in sorted(byELabel.keys()):
            if "good" in l: continue
            stanzas = byELabel[l]
            labelTracks.append("""
track {tn}_el
parent super_byLabel
compositeTrack on
shortLabel {shortL}
longLabel {longL}
type bigWig 9 +
visibility full
maxHeightPixels 32:12:8
autoScale on
dragAndDrop subTracks
hoverMetadata on
darkerLabels on
""".format(tn = "compos_" + sanitize(l),
           shortL=makeShortLabel(l),
           longL=makeLongLabel(l + " (" + str(len(stanzas)) + " exps)")))

        l = ExclusionLabels[0]
        stanzas = byELabel[l]
        labelTracks.append("""
track {tn}_el
parent super_byLabel
compositeTrack on
shortLabel {shortL}
longLabel {longL}
type bigWig 9 +
visibility full
maxHeightPixels 32:12:8
autoScale on
dragAndDrop subTracks
hoverMetadata on
darkerLabels on
""".format(tn = "compos_" + sanitize(l),
           shortL=makeShortLabel(l),
           longL=makeLongLabel(l + " (" + str(len(stanzas)) + " exps)")))
            
        for labelTrack in labelTracks:
            f.write(labelTrack)
            f.write('\n')
        for t, stanzas in byELabel.items():
            for stanza in sorted(stanzas):
                f.write(stanza)
                f.write('\n')

    def _doView(self, f):
        def out(tissue, exp, fileID):
            parent = "compos_" + sanitize(tissue)
            ret = ""
            for line in outputLines(self.data(tissue, exp, fileID, parent), 1):
                ret += line
            return ret

        tissueTracks = []
        for idx, tup in enumerate(self.exps):
            expID, exp, t, fileID = tup
            print(idx, "of", len(self.exps), t, '\t\t',
                  expID, exp.biosample_term_name, exp.age_display)
            tracks.append(out(t, exp, fileID))

        for t, stanzas in tracks.items():
            tissueTracks.append("""
track {tn}_view
parent super_byView
compositeTrack on
shortLabel {shortL}
longLabel {longL}
type bigWig 9 +
visibility full
maxHeightPixels 32:12:8
autoScale on
dragAndDrop subTracks
hoverMetadata on
darkerLabels on
""".format(tn = "compos_" + sanitize(t),
           shortL=makeShortLabel(t),
           longL=makeLongLabel(t + " (" + str(len(stanzas)) + " exps)")))

        for tissueTrack in tissueTracks:
            f.write(tissueTrack)
            f.write('\n')
        for t, stanzas in byTissue.items():
            for stanza in stanzas:
                f.write(stanza)
                f.write('\n')

    def bigBedWeng(self, f):
        name = "Weng rDHS Score"
        stanza = """
	track WengrDHS
	bigDataUrl http://users.wenglab.org/purcarom/dnase/Weng-rDHS-Score.bigBed
	visibility dense
	type bigBed 5
	spectrum on
        shortLabel {shortL}
	longLabel {longL}
	darkerLabels on
        priority 1
""".format(shortL=makeShortLabel("Weng"),
           longL=makeLongLabel(name))
        f.write(stanza)
        f.write('\n')        

    def bigBedStam(self, f):
        name = "Stam Master List No Overlap"
        stanza = """
	track StamDHS
	bigDataUrl http://users.wenglab.org/purcarom/dnase/Stam-Master-List-WM20180313-NoOverlap.bigBed
	visibility dense
	type bigBed 5
        spectrum on
	shortLabel {shortL}
	longLabel {longL}
	darkerLabels on
        priority 2
""".format(shortL=makeShortLabel("Stam"),
           longL=makeLongLabel(name))
        f.write(stanza)
        f.write('\n')        
                
    def run(self):
        self.loadExps()

        superTrackNames = {"byTissue": "DNase Signal byTissue",
                           "byLabel": "DNase Signal by exclusion label"
                           #"byView": "DNase Signal view"
        }
        superTracks = []
        for trackName, name in superTrackNames.items():
            superTracks.append("""
track super_{t}
superTrack on show
shortLabel {shortL}
longLabel {longL}
""".format(t = trackName,
           shortL=makeShortLabel(name),
           longL=makeLongLabel(name)))
        
        fnp = os.path.join(os.path.dirname(__file__), '../www/hg38/trackDb.txt')
        with open(fnp, 'w') as f:
            for superTrack in superTracks:
                f.write(superTrack)
                f.write('\n')
            self.load_exclusion_labels()
            self.load_newgroup()
            self._doByTissue(f)
            self._doByExclusionLabel(f)
            self.bigBedWeng(f)
            self.bigBedStam(f)
        printWroteNumLines(fnp)

ds = DNaseStam()
ds.run()
