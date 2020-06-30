#ifndef MTM_EX3_PARTB_MATRIX_H
#define MTM_EX3_PARTB_MATRIX_H

#include <string>
#include <iostream>
#include "Auxiliaries.h"

namespace mtm {
    template<class T>
    class Matrix {
    private:
        T **matrix;
        Dimensions dimensions;

        void verifyDimension(Dimensions dims);

    public:
        explicit Matrix(Dimensions dims, T initial = T());

        Matrix(T initial);

        Matrix(const Matrix &matrix1);

        ~Matrix();

        Matrix &operator=(const Matrix &matrix1);

        static Matrix Diagonal(int dim, T t, T initial = T());

        int height() const;

        int width() const;

        int size() const;

        Matrix transpose() const;

        Matrix operator-() const;

        Matrix operator-(const Matrix matrix1) const;

        Matrix &operator+=(const T t);

        //Matrix operator+(const int x) const;

        const T operator()(const int row, const int col) const;

        T &operator()(const int row, const int col);

        Matrix<bool> operator<(T n) const;

        Matrix<bool> operator>(T t) const;

        Matrix<bool> operator<=(T t) const;

        Matrix<bool> operator>=(T t) const;

        Matrix<bool> operator==(T t) const;

        Matrix<bool> operator!=(T t) const;

        template<class Function>
        Matrix apply(Function func) const;

        class iterator {
        private:
            Dimensions dims;
            int index;
            T **matrix;
        public:
            iterator(Dimensions dims, int i, T **matrix);

            ~iterator() = default;

            iterator(const iterator &it);

            iterator &operator=(const iterator &it);

            T &operator*() const;

            iterator operator++();

            iterator operator++(int);

            bool operator==(const iterator &it) const;

            bool operator!=(const iterator &it) const;

        };

        class const_iterator {
        private:
            Dimensions dims;
            int index;
            T **matrix;
        public:
            const_iterator(Dimensions dims, int i, T **matrix);

            ~const_iterator() = default;

            const_iterator(const const_iterator &it);

            const_iterator &operator=(const const_iterator &it);

            const T operator*() const;

            const_iterator operator++();

            const_iterator operator++(int);

            bool operator==(const const_iterator &it) const;

            bool operator!=(const const_iterator &it) const;

        };

        iterator begin();

        iterator end();

        const_iterator begin() const;

        const_iterator end() const;

        class AccessIllegalElement : public std::exception {
        public:
            std::string what() {
                std::string str = "Mtm matrix error: An attempt to access an illegal element";
                return str;
            }

        };

    class IllegalInitialization : public std::exception {
        public:
            IllegalInitialization() = default;
            std::string what() {
                std::string str = "Mtm matrix error: Illegal initialization values";
                return str;
            }

        };

        class DimensionMismatch : public std::exception {
        private:
            Dimensions dimensions1;
            Dimensions dimensions2;
        public:

            DimensionMismatch(Dimensions dims1, Dimensions dims2)
                    : dimensions1(dims1), dimensions2(dims2) {
            }

            std::string what() {
                int row1 = dimensions1.getRow();
                int col1 = dimensions1.getCol();
                int row2 = dimensions2.getRow();
                int col2 = dimensions2.getCol();

                std::string str = "Mtm matrix error: Dimension mismatch: ";
                str += "(" + std::to_string(row1) + "," + std::to_string(col1) + ")";
                str += " (" + std::to_string(row2) + "," + std::to_string(col2) + ")";
                return str;
            }

        };
        friend std::ostream &operator<<(std::ostream &os, const Matrix &matrix1){
            const_iterator begin = matrix1.begin();
            const_iterator end = matrix1.end();
            return printMatrix(os, begin, end, matrix1.width());
        }


        friend Matrix<T> operator+(const Matrix<T> &matrix1, const Matrix<T> &matrix2) {
            if (matrix1.size() == 1) {
                Matrix<T> new_mat = matrix2;
                const T t = matrix1(0, 0);
                for (int i = 0; i < matrix2.height(); ++i) {
                    for (int j = 0; j < matrix2.width(); ++j) {
                        new_mat(i, j) = t + new_mat(i, j);
                    }
                }
                return new_mat;
            } else if (matrix2.size() == 1) {
                Matrix<T> new_mat = matrix1;
                const T t = matrix2(0, 0);
                for (int i = 0; i < matrix1.height(); ++i) {
                    for (int j = 0; j < matrix1.width(); ++j) {
                        new_mat(i, j) += t;
                    }
                }
                return new_mat;
            } else {
                if (matrix1.width() != matrix2.width() || matrix1.height() != matrix2.height()) {
                    throw Matrix<T>::DimensionMismatch(matrix1.dimensions, matrix2.dimensions);
                }
                Matrix<T> new_mat = matrix1;
                for (int i = 0; i < matrix1.height(); ++i) {
                    for (int j = 0; j < matrix1.width(); ++j) {
                        new_mat(i, j) += matrix2(i, j);
                    }
                }
                return new_mat;
            }
        }
    };


//    template<class T>
//    Matrix<T> operator+(const Matrix<T> &matrix1, const Matrix<T> &matrix2)

    template<class T>
    bool any(const Matrix<T> &matrix) {
        for (int i = 0; i < matrix.height(); ++i) {
            for (int j = 0; j < matrix.width(); ++j) {
                if (static_cast<bool>(matrix(i, j))) { // todo check if it works
                    return true;
                }
            }
        }
        return false;
    }

    template<class T>
    bool all(const Matrix<T> &matrix) {
        for (int i = 0; i < matrix.height(); ++i) {
            for (int j = 0; j < matrix.width(); ++j) {
                if (!(static_cast<bool>(matrix(i, j)))) { // todo check if it works
                    return false;
                }
            }
        }
        return true;
    }

/** Matrix constructor
 * initialize the fields and allocate new matrix with the desired type
 * @tparam T - The type which the matrix is based on
 * @param dims - The dimension of the matrix has to be Dimension type.
 * @param initial - The same type of T.
 */
    template<class T>
    Matrix<T>::Matrix(Dimensions dims, T initial) : matrix(nullptr), dimensions(dims) {
        //verifyDimension(dims);

        if (this->dimensions.getCol() <= 0 || this->dimensions.getRow() <= 0) {
            throw Matrix::IllegalInitialization();
        }
        int row = dimensions.getRow();
        int col = dimensions.getCol();
        matrix = new T *[row];
        matrix[0] = new T[row * col];
        for (int i = 1; i < row; ++i)
            matrix[i] = matrix[0] + i * col;
        for (int i = 0; i < row; i++) {
            for (int j = 0; j < col; j++) {
                matrix[i][j] = initial;
            }
        }

    }

/**
 * constructor which create a (1,1) dimension matrix with the initial T.
 * @tparam T - The desired type.
 * @param initial - The initial T which will be the only element of the matrix.
 */
    template<class T>
    Matrix<T>::Matrix(T initial) : matrix(nullptr), dimensions(1, 1) {
        int row = 1;
        int col = 1;
        matrix = new T *[row];
        matrix[0] = new T[row * col];
        for (int i = 1; i < row; ++i)
            matrix[i] = matrix[0] + i * col;
        for (int i = 0; i < row; i++) {
            for (int j = 0; j < col; j++) {
                matrix[i][j] = initial;
            }
        }
    }

/**
 * copy constructor - Initial matrix with the same data of other matrix , allocate a new memory
 * for the new matrix.
 * @tparam T - The desired type.
 * @param matrix1 - The matrix which will initialized with (from the same type).
 */
    template<class T>
    Matrix<T>::Matrix(const Matrix &matrix1) : matrix(nullptr), dimensions(matrix1.dimensions) {
        int row = dimensions.getRow();
        int col = dimensions.getCol();
        matrix = new T *[row];
        if (row) {
            matrix[0] = new T[row * col];
            for (int i = 1; i < row; ++i)
                matrix[i] = matrix[0] + i * col;
        }
        for (int i = 0; i < row; i++) {
            for (int j = 0; j < col; j++) {
                matrix[i][j] = matrix1.matrix[i][j];
            }
        }
    }

/**
 * d'tor - Free al the memory allocated in the object.
 * @tparam T - The desired type.
 */
    template<class T>
    Matrix<T>::~Matrix() {
        int row = this->dimensions.getRow();
        if (row)
            delete[] matrix[0];
        delete[] matrix;
    }

/**
 * operator= - assigment operator create a new matrix equals to matrix 1 and replace it with 'this' matrix.
 * @tparam T - The desired type.
 * @param matrix1 - matrix typed T.
 * @return - Matrix which equals to matrix1.
 */
    template<class T>
    Matrix<T> &Matrix<T>::operator=(const Matrix &matrix1) {
        int row = matrix1.dimensions.getRow();
        int col = matrix1.dimensions.getCol();
        T **temp = new T *[row];
        temp[0] = new T[row * col];
        for (int i = 1; i < row; ++i) {
            temp[i] = temp[0] + i * col;
        }
        delete[] matrix;
        matrix = temp;
        for (int i = 0; i < row; i++) {
            for (int j = 0; j < col; j++) {
                matrix[i][j] = matrix1.matrix[i][j];
            }
        }
        this->dimensions = matrix1.dimensions;
        return *this;
    }

    template<class T>
    void Matrix<T>::verifyDimension(Dimensions dims) {
        if (dims.getCol() <= 0 || dims.getRow() <= 0) {
            throw IllegalInitialization();
        }
    }

    template<class T>
    Matrix<T> Matrix<T>::Diagonal(int dim, T t, T initial) {
        if (dim <= 0) {
            throw IllegalInitialization();
        }
        Dimensions squared_dims(dim, dim);
        Matrix<T> result(squared_dims, initial);
        for (int i = 0; i < dim; ++i) {
            result.matrix[i][i] = t;
        }
        return result;
    }

    template<class T>
    int Matrix<T>::height() const {
        int result = dimensions.getRow();
        return result;
    }

    template<class T>
    int Matrix<T>::width() const {
        int result = dimensions.getCol();
        return result;
    }

    template<class T>
    int Matrix<T>::size() const {
        int result = (this->height()) * (this->width());
        return result;
    }

    template<class T>
    Matrix<T> Matrix<T>::transpose() const {
        Dimensions transposed_dimensions(this->dimensions.getCol(), this->dimensions.getRow());// in purpose in order to make transposed matrix
        Matrix<T> transposed(transposed_dimensions);
        int col = transposed.dimensions.getCol();
        int row = transposed.dimensions.getRow();
        for (int i = 0; i < row; ++i) {
            for (int j = 0; j < col; ++j) {
                transposed.matrix[i][j] = this->matrix[j][i];
            }
        }
        return transposed;
    }

    template<class T>
    Matrix<T> Matrix<T>::operator-() const {
        Matrix<T> result(this->dimensions);
        int col = this->dimensions.getCol();
        int row = this->dimensions.getRow();
        for (int i = 0; i < row; ++i) {
            for (int j = 0; j < col; ++j) {
                result.matrix[i][j] = -this->matrix[i][j];
            }
        }
        return result;
    }

    template<class T>
    Matrix<T> Matrix<T>::operator-(const Matrix matrix1) const {
        if (this->dimensions != matrix1.dimensions) {
            throw Matrix<T>::DimensionMismatch(this->dimensions, matrix1.dimensions);
        }        Matrix<T> temp = -matrix1;
        Matrix result = *this + temp;
        return result;
    }

    template<class T>
    Matrix<T> &Matrix<T>::operator+=(const T t) {
        int row = this->dimensions.getRow();
        int col = this->dimensions.getCol();
        for (int i = 0; i < row; ++i) {
            for (int j = 0; j < col; ++j) {
                this->matrix[i][j] += t;
            }
        }
        return *this;
    }


    template<class T>
    const T Matrix<T>::operator()(const int row, const int col) const {
        if (col < 0 || col >= this->dimensions.getCol() || row < 0 || row >= this->dimensions.getRow()) {
            throw AccessIllegalElement();
        }
        const T result = this->matrix[row][col];
        return result;
    }

    template<class T>
    T &Matrix<T>::operator()(const int row, const int col) {
        if (col < 0 || col >= this->dimensions.getCol() || row < 0 || row >= this->dimensions.getRow()) {
            throw AccessIllegalElement();
        }
        return this->matrix[row][col];
    }

    template<class T>
    Matrix<bool> Matrix<T>::operator<(T n) const {
        Matrix<bool> result(this->dimensions, false);
        int row = this->dimensions.getRow();
        int col = this->dimensions.getCol();
        for (int i = 0; i < row; ++i) {
            for (int j = 0; j < col; ++j) {
                if (this->matrix[i][j] < n) {
                    result(i, j) = true;
                }
            }
        }
        return result;

    }

    template<class T>
    Matrix<bool> Matrix<T>::operator>(T t) const {
        Matrix<bool> result(this->dimensions, false);
        int row = this->dimensions.getRow();
        int col = this->dimensions.getCol();
        for (int i = 0; i < row; ++i) {
            for (int j = 0; j < col; ++j) {
                if (this->matrix[i][j] > t) {
                    result(i, j) = true;
                }
            }
        }
        return result;
    }

    template<class T>
    Matrix<bool> Matrix<T>::operator<=(T t) const {
        Matrix<bool> result(this->dimensions, false);
        int row = this->dimensions.getRow();
        int col = this->dimensions.getCol();
        for (int i = 0; i < row; ++i) {
            for (int j = 0; j < col; ++j) {
                if (this->matrix[i][j] <= t) {
                    result(i, j) = true;
                }
            }
        }
        return result;
    }

    template<class T>
    Matrix<bool> Matrix<T>::operator>=(T t) const {
        Matrix<bool> result(this->dimensions, false);
        int row = this->dimensions.getRow();
        int col = this->dimensions.getCol();
        for (int i = 0; i < row; ++i) {
            for (int j = 0; j < col; ++j) {
                if (this->matrix[i][j] >= t) {
                    result(i, j) = true;
                }
            }
        }
        return result;
    }

    template<class T>
    Matrix<bool> Matrix<T>::operator==(T t) const {
        Matrix<bool> result(this->dimensions, false);
        int row = dimensions.getRow();
        int col = dimensions.getCol();
        for (int i = 0; i < row; ++i) {
            for (int j = 0; j < col; ++j) {
                if (this->matrix[i][j] == t) {
                    result(i, j) = true;
                }
            }
        }
        return result;
    }

    template<class T>
    Matrix<bool> Matrix<T>::operator!=(T t) const {
        Matrix<bool> result(this->dimensions, false);
        int row = this->dimensions.getRow();
        int col = this->dimensions.getCol();
        for (int i = 0; i < row; ++i) {
            for (int j = 0; j < col; ++j) {
                if (this->matrix[i][j] != t) {
                    result(i, j) = true;
                }
            }
        }
        return result;
    }


    template<class T>
    Matrix<T>::iterator::iterator(mtm::Dimensions dims, int i, T **matrix) : dims(dims), index(i), matrix(matrix) {
    }

    template<class T>
    Matrix<T>::iterator::iterator(const Matrix::iterator &it) :  dims(it.dims), index(it.index), matrix(it.matrix) {
    }

    template<class T>
    typename Matrix<T>::iterator &Matrix<T>::iterator::operator=(const Matrix::iterator &it) {
        if (this == &it) {
            return *this;
        }
        dims = it.dims;
        matrix = it.matrix;
        index = it.index;
        return *this;
    }

    template<class T>
    T &Matrix<T>::iterator::operator*() const {
        if (index < 0 || index >= (dims.getCol() * dims.getRow())) {
            throw AccessIllegalElement();
        }
        int column = index % dims.getCol();
        int row = index / dims.getCol();
        return matrix[row][column];
    }

    template<class T>
    typename Matrix<T>::iterator Matrix<T>::iterator::operator++() {
        index++;
        return *this;
    }

    template<class T>
    typename Matrix<T>::iterator Matrix<T>::iterator::operator++(int) {
        iterator it = *this;
        this->index++;
        return it;
    }

    template<class T>
    bool Matrix<T>::iterator::operator==(const Matrix::iterator &it) const {
        return !operator!=(it);
    }

    template<class T>
    bool Matrix<T>::iterator::operator!=(const Matrix::iterator &it) const {
        bool b = ((dims != it.dims) || (index != it.index) || (matrix != it.matrix));
        return b;
    }


    ///const iterator implement:
    template<class T>
    Matrix<T>::const_iterator::const_iterator(Dimensions dims, int i, T **matrix)
            : dims(dims), index(i), matrix(matrix) {
    }

    template<class T>
    Matrix<T>::const_iterator::const_iterator(const Matrix::const_iterator &it)
            : dims(it.dims), index(it.index), matrix(it.matrix) {
    }

    template<class T>
    typename Matrix<T>::const_iterator &Matrix<T>::const_iterator::operator=(const Matrix::const_iterator &it) {
        if (this == &it) {
            return *this;
        }
        dims = it.dims;
        matrix = it.matrix;
        index = it.index;
        return *this;
    }

    template<class T>
    const T Matrix<T>::const_iterator::operator*() const {
        if (index < 0 || index >= (dims.getCol() * dims.getRow())) {
            throw AccessIllegalElement();
        }
        int column = index % dims.getCol();
        int row = index / dims.getCol();
        const T t = matrix[row][column];
        return t;
    }

    template<class T>
    typename Matrix<T>::const_iterator Matrix<T>::const_iterator::operator++() {
        this->index++;
        return *this;
    }

    template<class T>
    typename Matrix<T>::const_iterator Matrix<T>::const_iterator::operator++(int) {
        const_iterator it = *this;
        this->index++;
        return it;
    }

    template<class T>
    bool Matrix<T>::const_iterator::operator!=(const Matrix::const_iterator &it) const {
        bool b = ((dims != it.dims) || (index != it.index) || (matrix != it.matrix));
        return b;
    }

    template<class T>
    bool Matrix<T>::const_iterator::operator==(const Matrix::const_iterator &it) const {
        return !operator!=(it);
    }


    template<class T>
    typename Matrix<T>::iterator Matrix<T>::begin() {
        iterator it(this->dimensions, 0, this->matrix);
        return it;
    }

    template<class T>
    typename Matrix<T>::iterator Matrix<T>::end() {
        iterator it(this->dimensions, dimensions.getCol() * dimensions.getRow(), this->matrix);
        return it;
    }

    template<class T>
    typename Matrix<T>::const_iterator Matrix<T>::begin() const {
        const_iterator it(this->dimensions, 0, this->matrix);
        return it;
    }

    template<class T>
    typename Matrix<T>::const_iterator Matrix<T>::end() const {
        const_iterator it(this->dimensions, dimensions.getCol() * dimensions.getRow(), this->matrix);
        return it;
    }

    template<class T>
    template<class Function>
    Matrix<T> Matrix<T>::apply(Function func) const {
        Matrix<T> result = *this;
        int row = dimensions.getRow();
        int col = dimensions.getCol();
        for (int i = 0; i < row; ++i) {
            for (int j = 0; j < col; ++j) {
                result(i ,j) = func(result(i ,j));
            }
        }
        return result;
    }


}

#endif //MTM_EX3_PARTB_MATRIX_H
