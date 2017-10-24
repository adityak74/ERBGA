#include <stdio.h>
 
#define SIZE (70) /* amount of bits */
#define ARRAY_SIZE(x) (x/8+(!!(x%8)))
 
char get_bit(char *array, int index);
void toggle_bit(char *array, int index);

void toggle_bit(char *array, int index) {
    array[index / 8] ^= 1 << (index % 8);
}
 
char get_bit(char *array, int index) {
    return 1 & (array[index / 8] >> (index % 8));
}

int main(void) {
    /* initialize empty array with the right size */
    char x[ARRAY_SIZE(SIZE)][ARRAY_SIZE(SIZE)][ARRAY_SIZE(2)] = {0};
    fprintf(stderr, "%d\n", get_bit(x[0/8][0/8], 0) );
    toggle_bit(x[0/8][0/8], 0 );
    fprintf(stderr, "%d\n", get_bit(x[0/8][0/8] , 0) );

    fprintf(stderr, "%d\n", get_bit(x[0/8][0/8], 1) );
    toggle_bit(x[0/8][0/8], 1 );
    fprintf(stderr, "%d\n", get_bit(x[0/8][0/8] , 1) );

    fprintf(stderr, "ARRAY_SIZE : %d\n", ARRAY_SIZE(SIZE));
    fprintf(stderr, "Size : %lu\n", sizeof(x));
    fprintf(stderr, "Size : %d\n", ARRAY_SIZE(SIZE)*ARRAY_SIZE(SIZE)*ARRAY_SIZE(2));  

    return 0;
}