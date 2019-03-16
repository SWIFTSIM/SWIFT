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
    size_t value = *hashmap_lookup(&m, key);

    if(value != key) error("Incorrect value (%zu) found for key: %zu", value, key);
    //else message("Retrieved element, Key: %zu Value: %zu", key, value);
  }

  message("Checking for invalid key..."); 
  if(hashmap_lookup(&m, NUM_KEYS + 1) != NULL) error("Key: %d shouldn't exist or be created.", NUM_KEYS + 1);

  message("Freeing hash table...");
  hashmap_free(&m);

}
