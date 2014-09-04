/*This file is part of SciPAL.

    SciPAL is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    SciPAL is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with SciPAL.  If not, see <http://www.gnu.org/licenses/>.

Copyright  S. C. Kramer , J. Hagemann  2010 - 2014
*/


#ifndef cublas_BlancImporter_H
#define cublas_BlancImporter_H


#include <iomanip>
#include <iostream>
#include <cstdio>
#include <vector>
#include <cmath>

#include <QtGlobal>
#include <QFile>
#include <QTextStream>

#include <deal.II/base/subscriptor.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/vector.h>

#include <lac/cublas_BlancMatrix.h>
using namespace dealii;

class BlancImporter
{
public:
    BlancImporter(bool opt_diag = true);
    virtual ~BlancImporter();



    template<typename number> Vector<number> *               readVector( Vector<number> *vector, const int eltdim, const char *filename );
    template<typename number> BlancMatrix<number> * readMatrix(BlancMatrix<number> *matrix, SparsityPattern *pattern, const int eltdim, const char *filename);
    SparsityPattern *      readMatrixPattern(SparsityPattern *pattern, const char *filename);

   //! void setOptimizedDiagonal( bool opt_diag ){this->optimized_diagonal=opt_diag; }

private:
    FILE* fileOpen(const char* filename, const char* mode);
    unsigned long getRows(const char* filename, const char* calledFunction);
    //!    unsigned long NoOfRows;
    bool optimized_diagonal;

};

BlancImporter::BlancImporter(bool opt_diag)
    : optimized_diagonal(opt_diag)
{

}

BlancImporter::~BlancImporter()
{

}
template<typename number>
Vector<number>* BlancImporter::readVector(Vector<number> *vector, const int eltdim, const char *filename)
{
    int i = 0, j = 0;
    unsigned int eltnum = getRows(filename, __FUNCTION__);
    number value = 0.0;

    vector->reinit(eltnum*eltdim);


    FILE* istream = fileOpen(filename, "re");

    for ( i=0 ; i<eltnum ; i++ ) {
        for(j=0; j<eltdim; j++) {

        fscanf ( istream,
                 "%lf",
                 &value );



        (*vector)(i*eltdim+j) = value;

        }/*j*/
    } /*i*/

    vector->compress();


    return vector;
}

template<typename number>
BlancMatrix<number>* BlancImporter::readMatrix(BlancMatrix<number> *matrix, SparsityPattern *pattern, const int eltdim, const char *filename)
{

    int i = 0;
    int j = 0;
    int k = 0;
    int l = 0;

    std::vector<size_t> colidx_vector;
    int colidx = 0;
    number value = 0.0;

    FILE* istream = fileOpen(filename, "re");


    matrix->reinit(*pattern, eltdim);
    matrix->clear_inverse();

    for ( i=0 ; i<matrix->m() ; i++ ) {

        if ( fscanf ( istream,
                      "%*s" ) == EOF ) {
            return ( NULL );
        } /*if*/

        while ( getc(istream) != '\n' );

        colidx_vector.clear();


        //! damit ist sichergestellt, dass der richtige Eintrag am richtigen Ort in der Matrix landet
        for ( l=0 ; l<pattern->row_length(i) ; l++ ) {
            colidx_vector.push_back( pattern->column_number(i,l) );
        } /*j*/

        if(this->optimized_diagonal)  sort(colidx_vector.begin(), colidx_vector.end());


        for ( j=0 ; j<matrix->dimension() ; j++ ) {


            for ( l=0 ; l<pattern->row_length(i) ; l++ ) {
               // colidx = pattern->column_number(i,l);
               colidx = colidx_vector.at(l);

                for ( k=0 ; k<matrix->dimension() ; k++ ) {

                    if ( fscanf ( istream, "%lf", &value) != 1 ) {
                        return ( NULL );
                    } /*if*/


                    matrix->set(i, colidx,j,k, value);


                } /*k*/

            } /*l*/

        } /*j*/

    } /*i*/

//!     diag-opt. muss noch eingebaut werden
    //!        colidx_vector.clear();

    //!        //! ColIdx-Vektor aufbauen um die Indizes anschlieszend zu sortieren
    //!        //! damit ist sichergestellt, dass der richtige Eintrag am richtigen Ort in der Matrix landet
    //!        for ( j=0 ; j<pattern->row_length(i) ; ++j ) {
    //!            colidx_vector.push_back( pattern->column_number(i,j) );
    //!        } /*j*/

    //!        if(this->optimized_diagonal)  sort(colidx_vector.begin(), colidx_vector.end());

    //!        for ( j=0 ; j<pattern->row_length(i) ; ++j ) {
    //!            colidx = colidx_vector.at(j);

    //!            if ( fscanf ( istream,
    //!                          "%lf",
    //!                          &value ) != 1 ) {
    //!                return ( NULL );
    //!            } /*if*/


    //!            matrix->set(i,colidx, value);

    //!        } /*j*/
    //!    } /*i*/

    fclose(istream);


    return matrix;
}

SparsityPattern* BlancImporter::readMatrixPattern(SparsityPattern *pattern, const char *filename)
{

    unsigned long eltnum = getRows(filename, __FUNCTION__);



    size_t i = 0;
    size_t j = 0;
    unsigned long row = 0;
    unsigned long rowlen = 0;
    unsigned long colidx = 0;
    std::vector<unsigned int> rowlengths(eltnum);

    //! Datei öffnen um die Zeilenlängen zu ermitteln
    FILE* istream = fileOpen(filename,"re");

    if ( fscanf ( istream,
                  "row: rowlen: column indices:" ) == EOF ) {

        return ( NULL );
    } /*if*/


    for ( i=1 ; i<=eltnum ; i++ ) {

        if ( fscanf ( istream,
                      "%lu",
                      &row ) != 1 ) {
            std::cout << "bad format" << std::endl;

            return ( NULL );
        } /*if*/

        if ( row != i ) {
            std::cout << "invalid input: row " <<row<<", i " << i <<std::endl;
            return ( NULL );
        }/*if*/

        if ( fscanf ( istream,
                      "%lu",
                      &rowlen ) != 1 ) {
            std::cout << "invalid input: rowlen " <<rowlen <<std::endl;
            return ( NULL );
        } /*if*/

        rowlengths.at(i-1)=rowlen;


        for ( j=1 ; j<=rowlen ; j++ ) {

            if ( fscanf ( istream,
                          "%lu",
                          &colidx ) != 1 ) {
                std::cout << "invalid input: colidx " << colidx <<std::endl;

                return ( NULL );
            } /*if*/

        } /*j*/

    } /*i*/

    fclose(istream);

    i = 0;
    j = 0;
    row = 0;
    rowlen = 0;
    colidx = 0;

    //! SparsityPattern reinitialisieren
    pattern->reinit(eltnum,eltnum, rowlengths, this->optimized_diagonal);

    //! Wieder durch die Datei gehen um die Spaltenindizes rauszulesen.
    istream = fileOpen(filename,"re");

    if ( fscanf ( istream,
                  "row: rowlen: column indices:" ) == EOF ) {

        return ( NULL );
    } /*if*/

    for ( i=1 ; i<=eltnum ; i++ ) {

        if ( fscanf ( istream,
                      "%lu",
                      &row ) != 1 ) {

            return ( NULL );
        } /*if*/

        if ( row != i ) {

            return ( NULL );
        }/*if*/

        if ( fscanf ( istream,
                      "%lu",
                      &rowlen ) != 1 ) {

            return ( NULL );
        } /*if*/


        for ( j=1 ; j<=rowlen ; j++ ) {
            if ( fscanf ( istream,
                          "%lu",
                          &colidx ) != 1 ) {

                return ( NULL );
            } /*if*/
            pattern->add(i-1, colidx-1);


        } /*j*/

    } /*i*/
    fclose(istream);

    pattern->compress();

    return pattern;

}


FILE* BlancImporter::fileOpen(const char *filename, const char *mode)
{
    FILE* fstream      = NULL;
    fstream = fopen ( filename, mode);

    if ( fstream == NULL )
    {
        fprintf ( stderr,
                  "%sput file '%s' could not be opened\n",
                  (strcmp(mode,"re")? "Out": "In"),
                  filename );

//        perror ("The following error occurred");
//        printf( "Value of errno: %d\n", errno );
    }

    return fstream;
}


unsigned long BlancImporter::getRows(const char *filename, const char* calledFunction)
{
    //!    QFile file(filename);
    //!    unsigned long rows=0;

    //!    if(file.open(QIODevice::ReadOnly))
    //!    {
    //!        QTextStream stream(&file);
    //!        while(!stream.atEnd())
    //!            rows++;

    //!    }
    //!    file.close();


    //! Ermitteln wie viele Zeilen das Besetzungsmuster hat, in dem man die Anz. d. Zeilen der Datei zaehlt und anschlieszend 2 abzieht.
    system( QString("wc -l "+QString(filename) +"> rowCount.out").toStdString().c_str() );
    FILE* file = fopen("rowCount.out","re");
    unsigned long rows;
    fscanf ( file,"%lu",&rows );
    fclose(file);
    system("rm rowCount.out");

    if( QString(calledFunction).contains("Vector") )
        rows-=1;
    else if(QString(calledFunction).contains("Pattern"))
        rows-=2;




    //!qDebug() << "called Function: "<< calledFunction;
    //!qDebug() << "rows: " << rows;

    return rows;
}

#endif
