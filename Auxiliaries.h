#ifndef EX3_AUXILARIES_H
#define EX3_AUXILARIES_H
#include <string>
#include <iostream>
#include <ostream>


namespace mtm { 
    class Dimensions {
        int row, col;
    public:
        Dimensions( int row_t,  int col_t);
        std::string toString() const;
        bool operator==(const Dimensions& other) const;
        bool operator!=(const Dimensions& other) const;
        int getRow() const ;
        int getCol() const ;
    };

    std::string printMatrix(const int* matrix,const Dimensions& dim);

    template<class ITERATOR_T>
    std::ostream& printMatrix(std::ostream& os,ITERATOR_T begin,
                                ITERATOR_T end, unsigned int width){
        unsigned int row_counter=0;
        for (ITERATOR_T it= begin; it !=end; ++it) {
            if(row_counter==width){
                row_counter=0;
                os<< std::endl;
            }
            os <<*it<<" ";
            row_counter++;
        }
        os<< std::endl;
        return os;
    }
}



#endif //EX3_AUXILARIES_H
