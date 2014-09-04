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


#ifndef cublas_algorithms_H
#define cublas_algorithms_H

#include <iomanip>

#include <base/sign.h>

#include <lac/development/cublas_Vector.h>


namespace SciPAL {

    template<typename T>
    struct numerical_zero {
        static T Tol();
    };


    template<>
    inline float  numerical_zero<float>::Tol()  { return 1e-8; }

    template<>
    inline double numerical_zero<double>::Tol() { return 1e-16; }

    template<typename T>
    inline T numerical_zero<T>::Tol() {
        AssertThrow(false,
                    dealii::ExcMessage("Not implemented.") );
    }


    // @sect3{struct: algorithms}
    //!
    //! namespaces vertragen keine template-Argumente, also nehmen wir eine Struktur.
    template<typename T>
    struct algorithms {

        template<typename BW>
        static bool check_orthogonality (SciPAL::Matrix<T, BW> & Q)
        {
            typedef SciPAL::Matrix<T, BW> M;

            if(Q.n_rows() != Q.n_cols()) return false;

            M Q_Q_t = Q * transpose<M>(Q);

            M I = dealii::IdentityMatrix(Q.n_cols());

            //! Q_Q_t.print();

            Q_Q_t -= I;

            T l2_error = Q_Q_t.l2_norm();

            if (l2_error < Q.n_cols()*numerical_zero<T>::Tol() )
                return true;
            else {
                std::cout << "Deviation from unit matrix : " << std::setprecision(16) << l2_error << "\n" << std::endl;
                return false;
            }
        }



        template <typename VecView>
        static void test_householder_property(int t, T sigma, T alpha,
                                              const Vector<T, typename VecView::BW> & v,
                                              const VecView & sub_col_A);

    template <typename VecView>
        inline static void compute_householder_vector(int t, T & alpha, T & sigma, Vector<T, typename VecView::BW> & W,
                                                      const T A_tt, const VecView & entries2kill)
        {

        typedef typename VecView::BW BW;

    #ifdef HOUSEHOLDER_DEBUG
        std::cout << "\nZu eliminierende Eintraege :" << std::endl; entries2kill.print();
    #endif
                //! Zwischendurch wird immer wieder ein Vektor gebraucht,
                //! dessen Elemente bis auf das erste Null sind. Dafuer missbrauchen
                //! wir die FullMatrixAccessor-Klasse, um Zugriff auf
                //! das Rohdatenarray zu haben.
                //! Bei einem dealii::Vector haetten wir das nicht.
                //! Noch eine Wrapper-Klasse zu schreiben habe ich keine Lust.
            //! XXX //! FullMatrixAccessor<T> zero_vector(W.size(), 1);

                //! $ \alpha = sign (A_{tt})\|y\|_2$
            alpha    =  sign(A_tt) * entries2kill.l2_norm();
            T beta_t = A_tt + alpha;
            sigma    = + beta_t/alpha;

    #ifdef HOUSEHOLDER_DEBUG
        std::cout << "t : " << t << ", alpha : " << std::setprecision(5)
            << alpha <<", y + alpha : "<< beta_t <<", sigma : "<< sigma << std::endl << std::endl;
    #endif
/* XXX
            zero_vector(t,0) = alpha;
            Vector<T, BW> alpha_v_m_d(zero_vector, 0, 0);
            zero_vector(t,0) = 0.;
*/
    #ifdef HOUSEHOLDER_DEBUG
            std::cout << "\nalpha-Vektor :" << std::endl;  alpha_v_m_d.print();
    #endif

            W = entries2kill;

           //! W.add(t ???, alpha); //! W += alpha_v_m_d;
            W *= 1./beta_t;

                //! An dieser Stelle ist $\hat u^{(i)}$ fertig
    #ifdef HOUSEHOLDER_DEBUG
            std::cout << "\nAn dieser Stelle ist der Householder-Vektor fertig :" << std::endl;
            W.print();
    #endif

    #ifdef HOUSEHOLDER_DEBUG
                //! TEST : $(I - \sigma uu^T)y = -(0, \ldots,0, \alpha, 0, \ldots,0)^T$ wobei
                //! $y$ die aktuell bearbeitete Spalte von $A$ ist.
                //! Kann analog fuer Zeilen gemacht werden.
            test_householder_property(t, sigma, -alpha, W, entries2kill);
    #endif
        }
    };
}





// @sect4{Funktion: test_householder_property}
//!
template<typename T>
template <typename VecView>
void SciPAL::algorithms<T>::test_householder_property(int t, T sigma, T alpha,
                                                       const SciPAL::Vector<T, typename VecView::BW> &v,
                                                       const VecView &col_A)
{

    typedef typename VecView::BW BW;

    int n_entries = v.size();

    SciPAL::Vector<T, BW> entries2kill(n_entries);

    entries2kill = col_A;

#ifdef HOUSEHOLDER_DEBUG
    SciPAL::Matrix<T> UUt(n_entries, n_entries);

    UUt.add_scaled_outer_product(-sigma, v, v);

    std::cout << "\nAn dieser Stelle hat man -sigma vv^T :" << std::endl;
    UUt.print();
#endif

    SciPAL::Matrix<T, BW> tmp_Id(n_entries, n_entries);
    tmp_Id = dealii::IdentityMatrix(n_entries);

    tmp_Id.add_scaled_outer_product(-sigma, v, v);

    SciPAL::Vector<T, BW>       new_col(n_entries);

    tmp_Id.vmult(new_col, entries2kill);

#ifdef HOUSEHOLDER_DEBUG
    std::cout << "\nAn dieser Stelle hat man XXX -("
            << -alpha << ", 0, ... , 0) :" << std::endl;
    new_col.print();
#endif
}

#endif
