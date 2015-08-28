#include <base/print_info.h>

void print_expr_info(QString expr_name)
{
QStringList names = expr_name.split("[with");
std::cout<<"evaluating expression :"<<names.at(1).toStdString()<<"\n"<<std::endl;
}
