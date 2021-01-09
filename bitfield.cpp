//
//  bitfield.cpp
//  dummylib
//
//  Created by Anuj Jain on 11/19/20.
//

#include "bitfield.hpp"

bit_queue::bit_queue(int max_items)        //Constructor
{
    int num_bit_arrays, i;
    num_bit_arrays = (max_items / BIT_ARRAY_SIZE) + 1;
    bit_arrays = new unsigned long [num_bit_arrays];
    if (bit_arrays == NULL) {
        fprintf(stderr, "Out of memory!\n");
        //exit(0);
    }
    for (i = 0; i < num_bit_arrays; i++)
        bit_arrays[i] = 0;
    max_bit_arrays = num_bit_arrays;
    num_items = 0;
}

void bit_queue::create_bit_queue( int max_items)
{
    int num_bit_arrays, i;
    num_bit_arrays = (max_items / BIT_ARRAY_SIZE) + 1;
    bit_arrays = new unsigned long [num_bit_arrays];
    if (bit_arrays == NULL) {
        fprintf(stderr, "Out of memory!\n");
        //exit(0);
    }
    for (i = 0; i < num_bit_arrays; i++)
        bit_arrays[i] = 0;
    max_bit_arrays = num_bit_arrays;
    num_items = 0;
}

int bit_queue::empty_bit_queue()
{
    //int i;
//    memset( bit_arrays, 0, sizeof(unsigned long)* max_bit_arrays);
    for (int i = 0; i < max_bit_arrays; i++)
    {
        bit_arrays[i] = 0;
    }
    num_items = 0;
    return 0;
}

unsigned long bit_queue::check_bit_obj_present(int obj )
{
    int index_bit_array = obj/BIT_ARRAY_SIZE;       //index of the bit array to use.
    int set_bit = obj - (index_bit_array * BIT_ARRAY_SIZE);     //The index of bit to be set in the chosen bit array
    unsigned long number_bitset = 0x1 << set_bit;        //The number with the required bit set.
    if ( (set_bit > BIT_ARRAY_SIZE) || (index_bit_array >= max_bit_arrays))
    {
        fprintf(stderr, "Bit manip problem!\n");
//        exit(0);
    }
    //will return non-zero if the bit is already set else zero.
    return (bit_arrays[index_bit_array] & number_bitset);
}

int bit_queue::bit_queue_pop(int obj )
{
    int index_bit_array = obj/BIT_ARRAY_SIZE;       //index of the bit array to use.
    int set_bit = obj - (index_bit_array * BIT_ARRAY_SIZE);     //The index of bit to be set in the chosen bit array
    unsigned long number_bit_unset = 0x1 << set_bit;        //The number with the required bit set.
    if (!(bit_arrays[index_bit_array] & number_bit_unset))
        return 1;           //bit not in the queue so nothing to do.
    number_bit_unset = (~number_bit_unset) & MAX_DEB_SEQ_SIZE;
    if ( (set_bit > BIT_ARRAY_SIZE) || (index_bit_array >= max_bit_arrays) )
    {
        fprintf(stderr, "Bit manip problem!\n");
//        exit(0);
    }
    bit_arrays[index_bit_array] &= number_bit_unset;  //set the bit corresponding to obj as 0
    num_items--;
    return 1;
}

int bit_queue::queue_add_bit( int obj )
{
    int index_bit_array = obj/BIT_ARRAY_SIZE;       //index of the bit array to use.
    int set_bit = obj - (index_bit_array * BIT_ARRAY_SIZE);     //The index of bit to be set in the chosen bit array
    unsigned long number_bitset = 0x1 << set_bit;        //The number with the required bit set.
    if ( (set_bit > BIT_ARRAY_SIZE) || (index_bit_array >= max_bit_arrays) )
    {
        fprintf(stderr, "Bit manip problem!\n");
//        exit(0);
    }
    //Setting the required bit in the appropriate bit array.
    if (!(bit_arrays[index_bit_array] & number_bitset))     //Need to set if not present already.
    {
        bit_arrays[index_bit_array] |= number_bitset;
        num_items++;
    }
    return 1;
}

int bit_queue::numItems()
{
    return num_items;
}


void bit_queue::destroy_bit_queue()
{
    delete [] bit_arrays;
}

