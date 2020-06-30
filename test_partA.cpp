#include <functional>
#include <string>
#include <iostream>
#include "IntMatrix.h"

#define NUMBER_OF_TESTS 5

void test_print_matrix(std::string title,const mtm::IntMatrix& mat){
    std::cout << title<<":" << std::endl << mat;
}

bool test_1(){
    mtm::Dimensions dim(2,3);
    mtm::IntMatrix mat_1(dim);
    mtm::IntMatrix mat_2(dim,2);
    mtm::IntMatrix mat_3(mat_1);
    const mtm::IntMatrix mat_4 = mat_2;

    if (mat_1.size() != 6){
        return false;
    }

    if(mat_1.width()!=3){
        return false;
    }

    if(mat_4.height()!=2){
        return false;
    }
    return true;
}

bool test_2(){
    mtm::Dimensions dim(1 ,4);
    mtm::IntMatrix mat_1(dim,2);
    const mtm::IntMatrix mat_2 = mtm::IntMatrix::Identity(10);
    if (mat_2.size() != 100){
        return false;
    }
    const mtm::IntMatrix mat_1_T(mat_1.transpose());
    if(mat_1_T.width()!=1){
        return false;
    }
    const mtm::IntMatrix mat_3 = mat_1_T.transpose()+mat_1;
    const mtm::IntMatrix mat_4 = -mat_3;
    mtm::IntMatrix mat_5 = 6+mat_4-mat_1 +5;
    test_print_matrix("test_2, mat_5", mat_5+=1);
    return true;
}

bool test_3(){
    mtm::Dimensions dim(2 ,3);
    const mtm::IntMatrix mat_1 = mtm::IntMatrix::Identity(1);
    if(mat_1(0,0)!=1){
        return false;
    }
    mtm::IntMatrix mat_2(dim,1);
    int counter = 0;
    for(int i = 0;i<dim.getRow();++i){
        for(int j = 0;j<dim.getCol();++j){
            mat_2(i,j)=counter;
            ++counter;
        }
    }
    if(mat_2(1,2)!=5){
        return false;
    }
    return true;
}

bool test_4(){
    mtm::Dimensions dim(2 ,3);
    const mtm::IntMatrix mat_1 = mtm::IntMatrix::Identity(6);
    mtm::IntMatrix mat_2(dim,1);
    if(mtm::all(mat_1)){
        return false;
    }
    if(mtm::any(mat_2==0)){
        return false;
    }
    int counter = 0;
    for(int i = 0;i<dim.getRow();++i){
        for(int j = 0;j<dim.getCol();++j){
            mat_2(i,j)=counter;
            ++counter;
        }
    }
    return true;
}

bool test_5(){
    mtm::Dimensions dim(2 ,3);
    mtm::IntMatrix mat_1(dim);
    int counter = 0;
    for(mtm::IntMatrix::iterator it =mat_1.begin();it != mat_1.end();it++){
        *it=counter;
        ++counter;
    }
    const mtm::IntMatrix mat_2(mat_1);
    int sum_elements = 0;
    for(mtm::IntMatrix::const_iterator it =mat_2.begin();it != mat_2.end();++it){
        sum_elements+=*it;
    }
    if(sum_elements!=15){
        return false;
    }
    return true;
}

void run_test(std::function<bool()> test, std::string test_name){
    if(!test()){
        std::cout<<test_name<<" FAILED."<<std::endl;
        return;
    }
    std::cout<<test_name<<" SUCCEEDED."<<std::endl;

}

int main(){
    std::function<bool()> tests[NUMBER_OF_TESTS]= {test_1,test_2,test_3,test_4,test_5};
    for(int i=0;i<NUMBER_OF_TESTS;++i){
        run_test(tests[i],"Test "+std::to_string(i+1));
    }
}