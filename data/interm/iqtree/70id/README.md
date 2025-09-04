commands on 13K variant (70% seq sim - reassigned):

1. iqtree3 -s representatives_70id_reassign_alignemnt.phy -m LG -fast -T AUTO -n 10 -ninit 10
2. iqtree3 -s representatives_70id_reassign_alignemnt.phy -m LG -fast -T AUTO -alrt 1000
3. iqtree3 -s representatives_70id_reassign_alignemnt.phy -m LG -T AUTO -B 1000

commands on 1.2K variant (30% seq sim - aln with FAMSA 2):

1. iqtree3 -s aln_famsa2.phy -fast -T AUTO -alrt 1000
2. iqtree3 -s aln_famsa2.phy -m VT+R8 -T AUTO -B 1000
3. iqtree3 -s aln_famsa2.phy -m VT+R8 -T AUTO -B 1000 --nmax 10000 --undo
