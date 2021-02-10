/*
 * This file is part of the Micro Python project, http://micropython.org/
 *
 * The MIT License (MIT)
 *
 * Copyright (c) 2018 Scott Shawcroft for Adafruit Industries
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

#include "shared-bindings/displayio/Bitmap.h"

#include <string.h>

#include "py/runtime.h"

#include "math.h"

void common_hal_displayio_bitmap_construct(displayio_bitmap_t *self, uint32_t width,
    uint32_t height, uint32_t bits_per_value) {
    uint32_t row_width = width * bits_per_value;
    // align to size_t
    uint8_t align_bits = 8 * sizeof(size_t);
    if (row_width % align_bits != 0) {
        self->stride = (row_width / align_bits + 1);
    } else {
        self->stride = row_width / align_bits;
    }
    self->width = width;
    self->height = height;
    self->data = m_malloc(self->stride * height * sizeof(size_t), false);
    self->read_only = false;
    self->bits_per_value = bits_per_value;

    if (bits_per_value > 8 && bits_per_value != 16 && bits_per_value != 32) {
        mp_raise_NotImplementedError(translate("Invalid bits per value"));
    }

    // Division and modulus can be slow because it has to handle any integer. We know bits_per_value
    // is a power of two. We divide and mod by bits_per_value to compute the offset into the byte
    // array. So, we can the offset computation to simplify to a shift for division and mask for mod.

    self->x_shift = 0; // Used to divide the index by the number of pixels per word. Its used in a
                       // shift which effectively divides by 2 ** x_shift.
    uint32_t power_of_two = 1;
    while (power_of_two < align_bits / bits_per_value ) {
        self->x_shift++;
        power_of_two <<= 1;
    }
    self->x_mask = (1 << self->x_shift) - 1; // Used as a modulus on the x value
    self->bitmask = (1 << bits_per_value) - 1;

    self->dirty_area.x1 = 0;
    self->dirty_area.x2 = width;
    self->dirty_area.y1 = 0;
    self->dirty_area.y2 = height;
}

uint16_t common_hal_displayio_bitmap_get_height(displayio_bitmap_t *self) {
    return self->height;
}

uint16_t common_hal_displayio_bitmap_get_width(displayio_bitmap_t *self) {
    return self->width;
}

uint32_t common_hal_displayio_bitmap_get_bits_per_value(displayio_bitmap_t *self) {
    return self->bits_per_value;
}

uint32_t common_hal_displayio_bitmap_get_pixel(displayio_bitmap_t *self, int16_t x, int16_t y) {
    if (x >= self->width || x < 0 || y >= self->height || y < 0) {
        return 0;
    }
    int32_t row_start = y * self->stride;
    uint32_t bytes_per_value = self->bits_per_value / 8;
    if (bytes_per_value < 1) {
        size_t word = self->data[row_start + (x >> self->x_shift)];

        return (word >> (sizeof(size_t) * 8 - ((x & self->x_mask) + 1) * self->bits_per_value)) & self->bitmask;
    } else {
        size_t* row = self->data + row_start;
        if (bytes_per_value == 1) {
            return ((uint8_t*) row)[x];
        } else if (bytes_per_value == 2) {
            return ((uint16_t*) row)[x];
        } else if (bytes_per_value == 4) {
            return ((uint32_t*) row)[x];
        }
    }
    return 0;
}

void common_hal_displayio_bitmap_blit(displayio_bitmap_t *self, int16_t x, int16_t y, displayio_bitmap_t *source,
                            int16_t x1, int16_t y1, int16_t x2, int16_t y2, uint32_t skip_index, bool skip_index_none) {
    // Copy complete "source" bitmap into "self" bitmap at location x,y in the "self"
    // Add a boolean to determine if all values are copied, or only if non-zero
    // If skip_value is encountered in the source bitmap, it will not be copied.
    // If skip_value is `None`, then all pixels are copied.

    if (self->read_only) {
        mp_raise_RuntimeError(translate("Read-only object"));
    }

    // simplest version - use internal functions for get/set pixels
    for (int16_t i=0; i < (x2-x1) ; i++) {
        if ( (x+i >= 0) && (x+i < self->width) ) {
            for (int16_t j=0; j < (y2-y1) ; j++){
                if ((y+j >= 0) && (y+j < self->height) ) {
                    uint32_t value = common_hal_displayio_bitmap_get_pixel(source, x1+i, y1+j);
                    if ( (skip_index_none) || (value != skip_index) ) { // write if skip_value_none is True
                            common_hal_displayio_bitmap_set_pixel(self, x+i, y+j, value);
                    }
                }
            }
        }
    }
}

void common_hal_displayio_bitmap_blitfancy(displayio_bitmap_t *self, int16_t ox, int16_t oy,
                                            int16_t dest_clip0_x, int16_t dest_clip0_y,
                                            int16_t dest_clip1_x, int16_t dest_clip1_y,
                                            displayio_bitmap_t *source, int16_t px, int16_t py,
                                            int16_t source_clip0_x, int16_t source_clip0_y,
                                            int16_t source_clip1_x, int16_t source_clip1_y,
                                            float angle,
                                            float scale,
                                            uint32_t skip_index, bool skip_index_none) {

    // *self: destination bitmap
    // ox: the (ox, oy) destination point where the source (px,py) point is placed
    // oy:
    // dest_clip0: (x,y) is the corner of the clip window on the destination bitmap
    // dest_clip1: (x,y) is the other corner of the clip window of the destination bitmap
    // *source: the source bitmap
    // px: the (px, py) point of rotation of the source bitmap
    // py:
    // source_clip0: (x,y) is the corner of the clip window on the source bitmap
    // source_clip1: (x,y) is the other of the clip window on the source bitmap
    // angle: angle of rotation in radians, positive is clockwise
    // scale: scale factor
    // skip_index: color index that should be ignored (and not copied over)
    // skip_index_none: if skip_index_none is True, then all color indexes should be copied
    //                                                     (that is, no color indexes should be skipped)


    // Copy complete "source" bitmap into "self" bitmap at location x,y in the "self"
    // Add a boolean to determine if all values are copied, or only if non-zero
    // If skip_value is encountered in the source bitmap, it will not be copied.
    // If skip_value is `None`, then all pixels are copied.


    // # Credit from https://github.com/wernsey/bitmap
    // # MIT License from
    // #  * Copyright (c) 2017 Werner Stoop <wstoop@gmail.com>
    // #
    // # *
    // # * #### `void bm_rotate_blit(Bitmap *dst, int ox, int oy, Bitmap *src, int px, int py, double angle, double scale);`
    // # *
    // # * Rotates a source bitmap `src` around a pivot point `px,py` and blits it onto a destination bitmap `dst`.
    // # *
    // # * The bitmap is positioned such that the point `px,py` on the source is at the offset `ox,oy` on the destination.
    // # *
    // # * The `angle` is clockwise, in radians. The bitmap is also scaled by the factor `scale`.
    // #
    // # void bm_rotate_blit(Bitmap *dst, int ox, int oy, Bitmap *src, int px, int py, double angle, double scale);


    // #     /*
    // #    Reference:
    // #    "Fast Bitmap Rotation and Scaling" By Steven Mortimer, Dr Dobbs' Journal, July 01, 2001
    // #    http://www.drdobbs.com/architecture-and-design/fast-bitmap-rotation-and-scaling/184416337
    // #    See also http://www.efg2.com/Lab/ImageProcessing/RotateScanline.htm
    // #    */


    if (self->read_only) {
        mp_raise_RuntimeError(translate("Read-only object"));
    }

    int16_t x,y;

    int16_t minx = dest_clip1_x;
    int16_t miny = dest_clip1_y;
    int16_t maxx = dest_clip0_x;
    int16_t maxy = dest_clip0_y;

    float sinAngle = sinf(angle);
    float cosAngle = cosf(angle);

    float dx, dy;

    /* Compute the position of where each corner on the source bitmap
    will be on the destination to get a bounding box for scanning */
    dx = -cosAngle * px * scale + sinAngle * py * scale + ox;
    dy = -sinAngle * px * scale - cosAngle * py * scale + oy;
    if(dx < minx) minx = (int16_t)dx;
    if(dx > maxx) maxx = (int16_t)dx;
    if(dy < miny) miny = (int16_t)dy;
    if(dy > maxy) maxy = (int16_t)dy;

    dx = cosAngle * (source->width - px) * scale + sinAngle * py * scale + ox;
    dy = sinAngle * (source->width - px) * scale - cosAngle * py * scale + oy;
    if(dx < minx) minx = (int16_t)dx;
    if(dx > maxx) maxx = (int16_t)dx;
    if(dy < miny) miny = (int16_t)dy;
    if(dy > maxy) maxy = (int16_t)dy;

    dx = cosAngle * (source->width - px) * scale - sinAngle * (source->height - py) * scale + ox;
    dy = sinAngle * (source->width - px) * scale + cosAngle * (source->height - py) * scale + oy;
    if(dx < minx) minx = (int16_t)dx;
    if(dx > maxx) maxx = (int16_t)dx;
    if(dy < miny) miny = (int16_t)dy;
    if(dy > maxy) maxy = (int16_t)dy;

    dx = -cosAngle * px * scale - sinAngle * (source->height - py) * scale + ox;
    dy = -sinAngle * px * scale + cosAngle * (source->height - py) * scale + oy;
    if(dx < minx) minx = (int16_t)dx;
    if(dx > maxx) maxx = (int16_t)dx;
    if(dy < miny) miny = (int16_t)dy;
    if(dy > maxy) maxy = (int16_t)dy;

    /* Clipping */
    if(minx < dest_clip0_x) minx = dest_clip0_x;
    if(maxx > dest_clip1_x - 1) maxx = dest_clip1_x - 1;
    if(miny < dest_clip0_y) miny = dest_clip0_y;
    if(maxy > dest_clip1_y - 1) maxy = dest_clip1_y - 1;

    float dvCol = cosAngle / scale;
    float duCol = sinAngle / scale;

    float duRow = dvCol;
    float dvRow = -duCol;

    float startu = px - (ox * dvCol + oy * duCol);
    float startv = py - (ox * dvRow + oy * duRow);

    float rowu = startu + miny * duCol;
    float rowv = startv + miny * dvCol;

    for(y = miny; y <= maxy; y++) {
        float u = rowu + minx * duRow;
        float v = rowv + minx * dvRow;
        for(x = minx; x <= maxx; x++) {
            if(u >= source_clip0_x  && u < source_clip1_x && v >= source_clip0_y && v < source_clip1_y) {
                uint32_t c = common_hal_displayio_bitmap_get_pixel(source, u, v);
                if( (skip_index_none) || (c != skip_index) ) {
                    common_hal_displayio_bitmap_set_pixel(self, x, y, c);
                }
            }
            u += duRow;
            v += dvRow;
        }
        rowu += duCol;
        rowv += dvCol;
    }
}

void common_hal_displayio_bitmap_set_pixel(displayio_bitmap_t *self, int16_t x, int16_t y, uint32_t value) {
    if (self->read_only) {
        mp_raise_RuntimeError(translate("Read-only object"));
    }
    // Update the dirty area.
    if (self->dirty_area.x1 == self->dirty_area.x2) {
        self->dirty_area.x1 = x;
        self->dirty_area.x2 = x + 1;
        self->dirty_area.y1 = y;
        self->dirty_area.y2 = y + 1;
    } else {
        if (x < self->dirty_area.x1) {
            self->dirty_area.x1 = x;
        } else if (x >= self->dirty_area.x2) {
            self->dirty_area.x2 = x + 1;
        }
        if (y < self->dirty_area.y1) {
            self->dirty_area.y1 = y;
        } else if (y >= self->dirty_area.y2) {
            self->dirty_area.y2 = y + 1;
        }
    }

    // Update our data
    int32_t row_start = y * self->stride;
    uint32_t bytes_per_value = self->bits_per_value / 8;
    if (bytes_per_value < 1) {
        uint32_t bit_position = (sizeof(size_t) * 8 - ((x & self->x_mask) + 1) * self->bits_per_value);
        uint32_t index = row_start + (x >> self->x_shift);
        uint32_t word = self->data[index];
        word &= ~(self->bitmask << bit_position);
        word |= (value & self->bitmask) << bit_position;
        self->data[index] = word;
    } else {
        size_t* row = self->data + row_start;
        if (bytes_per_value == 1) {
            ((uint8_t*) row)[x] = value;
        } else if (bytes_per_value == 2) {
            ((uint16_t*) row)[x] = value;
        } else if (bytes_per_value == 4) {
            ((uint32_t*) row)[x] = value;
        }
    }
}

displayio_area_t* displayio_bitmap_get_refresh_areas(displayio_bitmap_t *self, displayio_area_t* tail) {
    if (self->dirty_area.x1 == self->dirty_area.x2) {
        return tail;
    }
    self->dirty_area.next = tail;
    return &self->dirty_area;
}

void displayio_bitmap_finish_refresh(displayio_bitmap_t *self) {
    self->dirty_area.x1 = 0;
    self->dirty_area.x2 = 0;
}

void common_hal_displayio_bitmap_fill(displayio_bitmap_t *self, uint32_t value) {
    if (self->read_only) {
        mp_raise_RuntimeError(translate("Read-only object"));
    }
    // Update the dirty area.
    self->dirty_area.x1 = 0;
    self->dirty_area.x2 = self->width;
    self->dirty_area.y1 = 0;
    self->dirty_area.y2 = self->height;

    // build the packed word
    uint32_t word = 0;
    for (uint8_t i=0; i<32 / self->bits_per_value; i++) {
        word |= (value & self->bitmask) << (32 - ((i+1)*self->bits_per_value));
    }
    // copy it in
    for (uint32_t i=0; i<self->stride * self->height; i++) {
        self->data[i] = word;
    }
}
