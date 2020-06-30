#include <iostream>
#include <string>
#include "Matrix.h"

class Square {
public:
    int operator()(int val){
        return val*val;
    }
};

int main(){
    mtm::Dimensions dim_1(2,3);
    mtm::Dimensions dim_2(-2,3);
    mtm::Dimensions dim_3(3,2);
    try{
        const mtm::Matrix<std::string> mat_1(dim_2,"str1");
    } catch(mtm::Matrix<std::string>::IllegalInitialization& e){
        std::cout<<e.what()<<std::endl;
    }
    try{
        mtm::Matrix<std::string> mat_1(dim_1,"str1");
        const mtm::Matrix<std::string> mat_2((mat_1+std::string("+"))+mat_1);
        const mtm::Matrix<std::string> mat_3(mat_1.transpose());
        std::cout<<mat_2;
        mtm::Matrix<std::string> mat_4=mat_1+mat_3;
    } catch(mtm::Matrix<std::string>::DimensionMismatch& e){
        std::cout<<e.what()<<std::endl;
    }
    try{
        const mtm::Matrix<int> mat_1 = mtm::Matrix<int>::Diagonal(2,1);
        mtm::Matrix<int> mat_2 = mtm::Matrix<int>::Diagonal(2,-1);
        std::cout<<mtm::any(mat_1)<<" "<<mtm::any(mat_1-mat_2)<<" "<<mtm::all(mat_1+mat_2+1)<<std::endl;
        std::cout<<(-mat_2)(2,2)<<std::endl;
    } catch(mtm::Matrix<int>::AccessIllegalElement& e){
        std::cout<<e.what()<<std::endl;
    }
    try{
        mtm::Matrix<int> mat_1(dim_1);
        int counter=0;
        mtm::Matrix<int>::iterator it=mat_1.begin();
        for(;it!=mat_1.end();++it){
            *it=counter++;
        }
        std::cout<<(mat_1>=3);
        it++;
        std::cout<<(mat_1==3);
        *it;
    } catch(mtm::Matrix<int>::AccessIllegalElement& e){
        std::cout<<e.what()<<std::endl;
    }
    try{
        const mtm::Matrix<int> mat_1(dim_1,4);
        const mtm::Matrix<int> mat_2=mat_1.apply(Square());
        std::cout<<mat_1;
        std::cout<<mat_2(1,2)<<std::endl;
        mtm::Matrix<int>::const_iterator it = mat_2.end();
        std::cout<<*it;
    } catch(mtm::Matrix<int>::AccessIllegalElement& e){
        std::cout<<e.what()<<std::endl;
    }
}
