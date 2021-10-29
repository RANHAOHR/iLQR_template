#ifndef MATRIXOPERATION_H
#define MATRIXOPERATION_H

#include <vector>
#include <cstdio>
#include <iostream>

#include <cstdlib>
#include <cstdint>
#include <cmath>

std::vector< std::vector<double> > mult_ABC ( std::vector< std::vector<double> > &in_A, std::vector< std::vector<double> > &in_B, std::vector< std::vector<double> > &in_C);
std::vector< std::vector<double> > add_matrix(std::vector< std::vector<double> >  &in_A, std::vector< std::vector<double> > &in_B);
std::vector< std::vector<double> > mult_AB ( std::vector< std::vector<double> > &in_A, std::vector< std::vector<double> > &in_B);
std::vector<double> diff_vector(std::vector< double >  &in_A, std::vector<double> &in_B);
std::vector<double> add_vector(std::vector< double >  &in_A, std::vector<double> &in_B);
std::vector<double> left_Mult_AB (std::vector< std::vector<double> >  &in_A, std::vector<double> &in_B);
std::vector<double> right_Mult_AB ( std::vector<double> &in_A, std::vector< std::vector<double> > &in_B);
double dot_product ( std::vector<double>  &in_A, std::vector<double> &in_B);
std::vector< std::vector<double> > transpose_matrix(std::vector< std::vector<double> > &in_A );
bool inverse(std::vector< std::vector<double> > &in_A, std::vector< std::vector<double> > &out_inverse);
void adjoint(std::vector< std::vector<double> > &in_A, std::vector< std::vector<double> > &adj);
double determinant(std::vector< std::vector<double> > &in_A, int n);
void getCofactor( std::vector< std::vector<double> > &in_A, std::vector< std::vector<double> > &out_, int p, int q, int n);
double norm(std::vector<double> &in_A);
double norm(std::vector<double> &in_A, std::vector<double> &in_B);
double norm(std::vector<std::vector<double>> &in_A);

#endif //MATRIXOPERATION_H
