//
//  bitfield.hpp
//  dummylib
//
//  Created by Anuj Jain on 11/19/20.
//

#ifndef bitfield_hpp
#define bitfield_hpp

#include <stdio.h>

#define BIT_ARRAY_SIZE 31
#define MAX_DEB_SEQ_SIZE 0xFFFFFFFF

class bit_queue {
private:
    unsigned long *bit_arrays;
    //    int maxitems;
    int max_bit_arrays;
    int num_items;

public:
    //bit queues
    bit_queue(int maxitems);        //Constructor
    void create_bit_queue( int maxitems );
    int queue_add_bit( int obj );
    //int queue_pop_bit( bit_queue *q, int *result );
    //int bit_queue_has_items(bit_queue *q);
    int bit_queue_pop( int obj );
    unsigned long check_bit_obj_present(int obj );
    int empty_bit_queue();
    void destroy_bit_queue();
    int numItems();
};

#endif /* bitfield_hpp */
