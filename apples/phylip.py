import tempfile

def tophylip(tags,seqs):
    assert tags and seqs
    phy_fp = tempfile.NamedTemporaryFile(delete=False, mode='w+t')
    phy_fp.write("%d %d\n" % (len(tags), len(seqs[0])))
    for tag,seq in zip(tags,seqs):
        phy_fp.write("%s %s\n" % (tag,seq))
    return phy_fp.name