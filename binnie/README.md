bridgebuilder - binnie
=============

[Sir Alexander Richardson Binnie] [1] was a Victorian-era civil engineer responsible for several major engineering projects, including the [Vauxhall Bridge] [2] over the River Thames in central London.

The binnie component of the BridgeBuilder system takes as input the original BAM (in old coordinate system order) and a bridge-aligned BAM (with reads in original-BAM order) with no unmapped reads and it outputs three BAMs:
   * BAM of unchanged reads (in original order) consisting of reads that do not map to the bridge
   * BAM of newly bridge-mapped reads (unordered) consisting of reads that were unmapped in the original BAM but have now been mapped to the bridge
   * BAM to-be-remapped (unordered and unmapped) consisting of all other reads



[1]: https://en.wikipedia.org/wiki/Alexander_Binnie      "Sir Alexander Richardson Binnie"
[2]: https://en.wikipedia.org/wiki/Vauxhall_Bridge       "Vauxhall Bridge"

