/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2020 Camila Correa (camila.correa@uva.nl)
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
#ifndef SWIFT_DARK_MATTER_IO_H
#define SWIFT_DARK_MATTER_IO_H

/**
 * @brief Specifies which dm particle fields to write to a dataset
 *
 * @param gparts The dm particle array.
 * @param list The list of i/o properties to write.
 * @param num_fields The number of i/o fields to write.
 */
INLINE static int sidm_write_gparts(const struct dmpart* dmparts,
                                              struct io_props* list) {
    
    /* List what we want to write */
    list[0] = io_make_output_field("SIDMevents", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f,
                                   dmparts, sidm_data.num_sidm, "Number of DM-DM collisions the particle has had");
    return 1;
    
}


#endif /* SWIFT_DARK_MATTER_IO_H */
