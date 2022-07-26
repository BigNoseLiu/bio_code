from autopvs1 import AutoPVS1
demo = AutoPVS1('13-113803407-G-A', 'hg19')
#demo2 = AutoPVS1('13-113149093-G-A', 'hg38')
if demo.islof:
	print(demo.hgvs_c, demo.hgvs_p, demo.consequence, demo.pvs1.criterion, demo.pvs1.strength_raw, demo.pvs1.strength)
