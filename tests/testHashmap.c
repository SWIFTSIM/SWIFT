/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2019 James Willis (james.s.willis@durham.ac.uk).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

/* Config parameters. */
#include "../config.h"

/* Local headers. */
#include "swift.h"

#include "../src/c_hashmap/hashmap.h"

#define NUM_KEYS (26 * 1000 * 1000)

int main(int argc, char *argv[]) {

  hashmap_t m;

  message("Initialising hash table...");
  hashmap_init(&m);
  
  message("Populating hash table...");
  for(hashmap_key_t key=0; key<NUM_KEYS; key++) {
    hashmap_value_t value;
    value.value_st = key;
    hashmap_put(&m, key, value);
  }

  message("Dumping hashmap stats.");
  hashmap_print_stats(&m);

  message("Retrieving elements from the hash table...");
  for(hashmap_key_t key=0; key<NUM_KEYS; key++) {
    hashmap_value_t value = *hashmap_lookup(&m, key);

    if(value.value_st != key) error("Incorrect value (%zu) found for key: %zu", value.value_st, key);
    //else message("Retrieved element, Key: %zu Value: %zu", key, value);
  }

  message("Checking for invalid key..."); 
  if(hashmap_lookup(&m, NUM_KEYS + 1) != NULL) error("Key: %d shouldn't exist or be created.", NUM_KEYS + 1);

  message("Checking hash table size..."); 
  if(m.size != NUM_KEYS) error("The no. of elements stored in the hash table are not equal to the no. of keys. No. of elements: %zu, no. of keys: %d", m.size, NUM_KEYS);
  
  message("Freeing hash table...");
  hashmap_free(&m);

}
