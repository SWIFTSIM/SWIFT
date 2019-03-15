
/* Local headers. */
#include "swift.h"

#include "../src/c_hashmap/hashmap.h"

#define INITIAL_NUM_CHUNKS 16
#define NUM_KEYS 50000

int main(int argc, char *argv[]) {

  hashmap_t m;

  message("Initialising hash table...");
  hashmap_init(&m);
  
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

  message("Checking for invalid key..."); 
  if(hashmap_get(&m, NUM_KEYS + 1, 0) != NULL) error("Key: %d shouldn't exist or be created.", NUM_KEYS + 1);

  message("Freeing hash table...");
  hashmap_free(&m);

}
