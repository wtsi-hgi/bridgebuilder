bridgebuilder
=============

BridgeBuilder efficiently remaps BAM/SAM reads to a new reference by first building a "bridge" reference, first mapping to that bridge, and then remapping only a subset of reads to the full new reference. 


The [BridgeBuilder System] [1] consists of several components, including [Baker] [2], [Binnie] [3], & [Brunel] [4] and also relies on [Samtools] [5] for SAM/BAM manipulations and is currently tested using [bwa] [6] for mapping (although it could potentially use other aligners as well).


[1]: docs/BridgeBuilderSystemDiagram.png           "BridgeBuilder System Diagram"
[2]: baker/README.md                               "Baker"
[3]: binnie/README.md                              "Binnie"
[4]: brunel/README.md                              "Brunel"
[5]: https://github.com/samtools/samtools          "Samtools"
[6]: https://github.com/lh3/bwa                    "Burrows-Wheeler Aligner"

