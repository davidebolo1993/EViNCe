import pybedtools
import sys
import os

if __name__ == '__main__':

	bedtruth=pybedtools.BedTool(os.path.abspath(sys.argv[1]))

	B1=''
	B2=''

	for el in bedtruth.sort():

		key1=el.chrom
		start1=str(el.start-500)
		end1=str(el.end+500)
		
		key2=el[3]
		start2=str(int(el[4])-500)
		end2=str(int(el[5])+500)

		if key1 != key2:

			B1+=key1+'\t'+start1+'\t'+end1+'\n'
			B2+=key2+'\t'+start2+'\t'+end2+'\n'


	B1b=pybedtools.BedTool(B1,from_string=True).as_intervalfile()
	B2b=pybedtools.BedTool(B2,from_string=True).as_intervalfile()

	bedtest=pybedtools.BedTool(os.path.abspath(sys.argv[2]))

	TP=0
	FP=0
	FN=0

	if len(bedtest) != 0:

		seen1=list()
		seen2=list()

		for el in bedtest.sort():

			key1=el.chrom
			start1=el.start-500
			end1=el.end+500
			q1=pybedtools.Interval(key1,start1,end1)

			key2=el[3]
			start2=int(el[4])
			end2=int(el[5])
			q2=pybedtools.Interval(key2,start2,end2)

			r1=B1b.all_hits(q1)
			r2=B1b.all_hits(q2)

			if len(r1) != 0 or len(r2) != 0:


				if len(r1) != 0:

					if r1[0] not in seen1:

						TP+=1
						seen1.append(r1[0])

				else:

					if r2[0] not in seen2:

						TP+=1
						seen2.append(r2[0])

			else:

				r1= B2b.all_hits(q1)
				r2=	B2b.all_hits(q2)

				if len(r1) != 0 or len(r2) != 0:

					if len(r1) != 0:

						if r1[0] not in seen1:

							TP+=1
							seen1.append(r1[0])

					else:

						if r2[0] not in seen2:

							TP+=1
							seen2.append(r2[0])

				else:

					FP+=1

		FN=len(pybedtools.BedTool(B1,from_string=True))-TP


		R=TP/(TP+FN)
		P=TP/(TP+FP)
		F1=2*(P*R)/(P+R)
		print(str(TP)+ '\t' +str(FP) +'\t' + str(FN) +'\t' + str(P) + '\t' + str(R) + '\t' + str(F1))

	else:

		print("NaN\tNaN\tNan\tNaN\tNaN\tNaN")


