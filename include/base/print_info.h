#ifndef PRINT_INFO_H
#define PRINT_INFO_H
#include <QString>
#include <QStringList>
#include <iostream>
//! function to print information about expression evaluation
__host__
void print_expr_info(QString expr_name)
{
QStringList names = expr_name.split("[with");
std::cout<<"evaluating expression :"<<names.at(1).toStdString()<<"\n"<<std::endl;
}

#endif // PRINT_INFO_H

