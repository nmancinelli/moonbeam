import pandas as pd

def main():

	filename = input('File name:' )

	fout = open(filename + '.txt' , 'w')

	names = ["Radius","Density","Vp","Vs","_","__"]
	df = pd.read_csv(filename + '.out', delim_whitespace = True, skiprows = 5, names = names)

	df.Radius = df.Radius * 1000.
	df.Density = df.Density * 1000.
	df.Vp = df.Vp * 1000.
	df.Vs = df.Vs * 1000.
	
	dffld = df.query('Vs == 0.0')
	noc = len(dffld)
	
	#because of added layer
	if noc > 0:
		noc=noc+1
	
	fout.write('#\n')
	fout.write('0 0 1 1\n')
	fout.write('%d %d %d\n' % (len(df)+1, 0, noc) )
	
	for i, row in df.iterrows():
		
		s='%10.2f %8.2f %8.2f %8.2f 0 0 0 0 0\n'%(row.Radius, row.Density, row.Vp, row.Vs)
		if row.Radius > 330000. and row.Radius < 335000:
			s2='%10.2f %8.2f %8.2f %8.2f 0 0 0 0 0\n'%(row.Radius, dlast, vplast, vslast)
			fout.write(s2)
		
		fout.write(s)
		dlast = row.Density
		vplast = row.Vp
		vslast = row.Vs
	
	fout.close()


main()
