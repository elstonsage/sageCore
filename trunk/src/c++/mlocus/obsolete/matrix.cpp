//=============================================================================
//  File:       "matrix.cpp"
//
//  Author:     Bob Steagall
//
//  History:    Version 1.00.00
//
//  Notes:    
//
//  Copyright (c) 1999 R.C. Elston
//  All Rights Reserved
//============================================================================
//
#ifdef _MSC_VER
    #include <app/SAGEconfig.h>
    #pragma hdrstop
    #pragma warning (disable : 4661)
#endif

#include <memory.h>
#include <assert.h>
#include <mlocus.new/matrix.h>

namespace SAGE   {
namespace MLOCUS {

template class sparse_matrix<double>;

void f()
{
    sparse_matrix<double>   smd;
}


//-----------------------------------------------------------------------------
//  Function:   matrix_base()
//
//  Purpose:    The sole constructor for the matrix_base class.
//-----------------------------------------------------------------------------
//
matrix_base::matrix_base()
{
    init();
}


//-----------------------------------------------------------------------------
//  Function:   ~matrix_base()
//-----------------------------------------------------------------------------
//
matrix_base::~matrix_base()
{
    destroy();
}


void
matrix_base::swap(matrix_base& M)
{
    std::swap(my_rows, M.my_rows);
    std::swap(my_cols, M.my_cols);
    std::swap(my_ctldata, M.my_ctldata);
    std::swap(my_rawdata, M.my_rawdata);
}


//-----------------------------------------------------------------------------
//  Function:   init()
//
//  Purpose:    This function initializes the internal variables for this class.
//-----------------------------------------------------------------------------
//
void
matrix_base::init()
{
    my_rows = my_cols = 0;
    my_ctldata = my_rawdata = 0;
}


//-----------------------------------------------------------------------------
//  Function:   create()
//
//  Purpose:    This function creates new instances of matrix_base objects.
//-----------------------------------------------------------------------------
//
void*
matrix_base::create(uint usize, uint rows, uint cols)
{
    uint    j;
    char*   V;
    char**  M;

    assert(rows > 0);
    assert(cols > 0);

    //- Allocate the array that will hold the data.
    //
    V = new char[usize*rows*cols];
    assert(V != 0);
    my_rawdata = (void*) V;
 
    //- Initialize our baby array.
    //
    memset(V, 0, usize*rows*cols);

    //- Allocate the array of pointers to the data that implement the matrix.
    //
    M = new char* [rows];
    assert(M != 0);
    my_ctldata = (void*) M;

    //- Set size fields.
    //
    my_rows = rows;
    my_cols = cols;

    //- Initialize the array of pointers to the data.
    //
    for (j = 0;  j < my_rows;  ++j)
    {
        M[j] = V + (j * usize * cols);
    }

    return (void*) M;
}


//-----------------------------------------------------------------------------
//  Function:   destroy()
//
//  Purpose:    This function deletes any allocated memory and zeros the 
//              correspondint pointer(s).
//-----------------------------------------------------------------------------
//
void
matrix_base::destroy()
{
    if (my_rawdata)
    {
        delete [] (char*) my_rawdata;
        my_rawdata = 0;
    }
    
    if (my_ctldata)
    {
        delete [] (char*) my_ctldata;
        my_ctldata = 0;
    }
}


//-----------------------------------------------------------------------------
//  Functions:  copy()
//
//  Purpose:    This function copys the raw bits of one matrix into this one.
//-----------------------------------------------------------------------------
//
void*
matrix_base::copy(uint usize, const matrix_base& MB)
{
    //- If the new size of the matrix is different from the current size,
    //  or data does not currently exist in this object, adjustments
    //  need to be made.
    //
    if (my_rows != MB.my_rows  ||  my_cols != MB.my_cols  ||  !my_rawdata)
    {
        destroy();
        create(usize, MB.my_rows, MB.my_cols);
    }
    
    //- Copy the raw bits into this object.
    //
    memcpy((char*) my_rawdata, (char*) MB.my_rawdata, usize*my_rows*my_cols);

    //- Return the matrix pointer taking into account the lower bound in the
    //  j-direction.
    //
    return (void*) my_ctldata;
}


//-----------------------------------------------------------------------------
//  Function:   matrix_base::modify()
//
//  Purpose:    This function modifies the size of this matrix.  If necessary, 
//              raw bits are copied from the old version of 'this'.
//-----------------------------------------------------------------------------
//
void*
matrix_base::modify(uint usize, uint new_rows, uint new_cols)
{
    assert(new_rows > 0  &&  new_cols > 0);

    //- Only change if a real change is requested.
    //
    if (new_rows != my_rows  ||  new_cols != my_cols)
    {
        //- Create a new matrix.  Afterward, we're going to cheat and
        //  pull values out directly of the new matrix to replace values
        //  in this matrix.
        //
        matrix_base   UBM;        

        UBM.create(usize, new_rows, new_cols);

        //- Copy the old data (if there was any) into the new matrix.  Delete
        //  the old data array.
        //
        if (my_rawdata)
        {
            //- Calculate the rowsize of the old and new matrices in bytes.
            //  Calculate the byte offset of each row.
            //
            uint    old_rowsize = usize*my_cols;
            uint    new_rowsize = usize*new_cols;

            //- Determine which is smaller: the old rowsize in bytes, or the
            //  new rowsize in bytes.
            //
            uint    rowsize = std::min(old_rowsize, new_rowsize);

            //- Determine the smaller of the old number of rows, or the new
            //  number of rows.
            //
            uint    colsize = std::min(my_cols, new_cols);

            //- Copy the old data into the new array, row by row.
            //
            uint new_offset = 0;
            uint old_offset = 0;

            for (uint j = 0; j < colsize; j++)
            {
                memcpy((char*) UBM.my_rawdata + new_offset,
                       (char*) my_rawdata + old_offset,
                       rowsize);
                
                new_offset += new_rowsize;
                old_offset += old_rowsize;
            }
        }

        //- Delete existing old control data and zero the internal pointers.
        //        
        destroy();

        //- Copy the values of the control data pointers from the temporary
        //  matrix into the control data pointers of this.
        //
        my_rawdata = UBM.my_rawdata;
        my_ctldata = UBM.my_ctldata;

        //- Set the new row and column sizes.
        //
        my_rows = new_rows;
        my_cols = new_cols;
        
        //- Re-initialize the values (specifically the data pointers) in the 
        //  temp matrix so that our new data is not deleted when the temp 
        //  matrix goes out of scope.
        //
        UBM.init();
    }

    //- Return the matrix pointer taking into account the lower bound in the
    //  j-direction.
    //
    return (void*) my_ctldata;
}

} // End namespace MLOCUS
} // End namespace SAGE

