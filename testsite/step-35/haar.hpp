/** The 1D Haar Transform **/
const double ONEOVERSQRTTWO=1.0/sqrt(2.0);
void haar1d(double *vec, int n)
{
	int i=0;
	int w=n;
	double *vecp = new double[n];
	for(i=0;i<n;i++)
		vecp[i] = 0;

	while(w>1)
	{
		w/=2;
		for(i=0;i<w;i++)
		{
			vecp[i] = (vec[2*i] + vec[2*i+1])*ONEOVERSQRTTWO;
			vecp[i+w] = (vec[2*i] - vec[2*i+1])*ONEOVERSQRTTWO;
		}

		for(i=0;i<(w*2);i++)
			vec[i] = vecp[i]; 
	}

	delete [] vecp;
}


/** A Modified version of 1D Haar Transform, used by the 2D Haar Transform function **/
void haar1(double *vec, int n, int w)
{
	int i=0;
	double *vecp = new double[n];
	for(i=0;i<n;i++)
		vecp[i] = 0;

	w/=2;
	for(i=0;i<w;i++)
	{
		vecp[i] = (vec[2*i] + vec[2*i+1])*ONEOVERSQRTTWO;
		vecp[i+w] = (vec[2*i] - vec[2*i+1])*ONEOVERSQRTTWO;
	}

	for(i=0;i<(w*2);i++)
		vec[i] = vecp[i];

	delete [] vecp;
}
void ihaar1(double *vec, int n, int w)
{
	int i=0;
	double *vecp = new double[n];
	for(i=0;i<n;i++)
		vecp[i] = 0;

	w/=2;
	for(i=0;i<w;i++)
	{
		vecp[2*i]=(vec[i]+vec[i+w])*ONEOVERSQRTTWO;
		vecp[2*i+1]=(vec[i]-vec[i+w])*ONEOVERSQRTTWO;
	}

	for(i=0;i<(w*2);i++)
		vec[i] = vecp[i];

	delete [] vecp;
}

/** The 2D Haar Transform **/
void haar2(double *matrix, int rows, int cols)
{
	double *temp_row = new double[cols];
	double *temp_col = new double[rows];

	int i=0,j=0;
	int w = cols, h=rows;
	while(w>1 || h>1)
	{
		if(w>1)
		{
			for(i=0;i<h;i++)
			{
				for(j=0;j<cols;j++)
					temp_row[j] = matrix[i*cols+j];

				haar1(temp_row,cols,w);

				for(j=0;j<cols;j++)
					matrix[i*cols+j] = temp_row[j];
			}
		}

		if(h>1)
		{
			for(i=0;i<w;i++)
			{
				for(j=0;j<rows;j++)
					temp_col[j] = matrix[j*cols+i];
				haar1(temp_col, rows, h);
				for(j=0;j<rows;j++)
					matrix[j*cols+i] = temp_col[j];
			}
		}

		if(w>1)
			w/=2;
		if(h>1)
			h/=2;
	}

	delete [] temp_row;
	delete [] temp_col;
}
void ihaar2(double *matrix, int rows, int cols)
{
	double *temp_row = new double[cols];
	double *temp_col = new double[rows];

	int i=0,j=0;
	int w = 2, h=2;
	while(w<=cols  || h <= rows)
	{
		if(w<=cols)
		{
			for(i=0;i<h;i++)
			{
				for(j=0;j<cols;j++)
					temp_row[j] = matrix[i*cols+j];

				ihaar1(temp_row,cols,w);

				for(j=0;j<cols;j++)
					matrix[i*cols+j] = temp_row[j];
			}
		}
		if(h<=rows)
		{
			for(i=0;i<w;i++)
			{
				for(j=0;j<rows;j++)
					temp_col[j] = matrix[j*cols+i];
				ihaar1(temp_col, rows, h);
				for(j=0;j<rows;j++)
					matrix[j*cols+i] = temp_col[j];
			}
		}

		if(w<=cols)
			w*=2;
		if(h<=rows)
			h*=2;
	}

	delete [] temp_row;
	delete [] temp_col;
}
