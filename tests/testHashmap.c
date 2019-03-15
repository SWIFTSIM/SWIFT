
/* Local headers. */
#include "swift.h"

#include "../src/c_hashmap/hashmap.h"

#define INITIAL_NUM_CHUNKS 16
#define NUM_KEYS 50000

int main(int argc, char *argv[]) {

  hashmap_t m;

  message("Initialising hash table...");
  hashmap_init(&m);
  
  message("Allocating chunks for the hash table...");
  hashmap_allocate_chunks(&m, INITIAL_NUM_CHUNKS);

  message("Populating hash table...");
  for(size_t key=0; key<NUM_KEYS; key++) {
    hashmap_put(&m, key, key);
  }

  message("Retrieving elements from the hash table...");
  for(size_t key=0; key<NUM_KEYS; key++) {
    size_t value = *hashmap_get(&m, key, 0);

    if(value != key) error("Incorrect value (%zu) found for key: %zu", value, key);
    //else message("Retrieved element, Key: %zu Value: %zu", key, value);
  }

  message("Freeing hash table...");
  hashmap_free(&m);

}
