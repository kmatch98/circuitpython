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

#include <stdint.h>

#include "lib/utils/context_manager_helpers.h"
#include "py/binary.h"
#include "py/objproperty.h"
#include "py/runtime.h"
#include "shared-bindings/microcontroller/Pin.h"
#include "shared-bindings/util.h"
#include "supervisor/shared/translate.h"

//| class Bitmap:
//|     """Stores values of a certain size in a 2D array"""
//|
//|     def __init__(self, width: int, height: int, value_count: int) -> None:
//|         """Create a Bitmap object with the given fixed size. Each pixel stores a value that is used to
//|         index into a corresponding palette. This enables differently colored sprites to share the
//|         underlying Bitmap. value_count is used to minimize the memory used to store the Bitmap.
//|
//|         :param int width: The number of values wide
//|         :param int height: The number of values high
//|         :param int value_count: The number of possible pixel values."""
//|         ...
//|
STATIC mp_obj_t displayio_bitmap_make_new(const mp_obj_type_t *type, size_t n_args, const mp_obj_t *pos_args, mp_map_t *kw_args) {
    mp_arg_check_num(n_args, kw_args, 3, 3, false);
    uint32_t width = mp_obj_get_int(pos_args[0]);
    uint32_t height = mp_obj_get_int(pos_args[1]);
    uint32_t value_count = mp_obj_get_int(pos_args[2]);
    uint32_t bits = 1;

    if (value_count == 0) {
        mp_raise_ValueError(translate("value_count must be > 0"));
    }
    while ((value_count - 1) >> bits) {
        if (bits < 8) {
            bits <<= 1;
        } else {
            bits += 8;
        }
    }

    displayio_bitmap_t *self = m_new_obj(displayio_bitmap_t);
    self->base.type = &displayio_bitmap_type;
    common_hal_displayio_bitmap_construct(self, width, height, bits);

    return MP_OBJ_FROM_PTR(self);
}
//|     width: int
//|     """Width of the bitmap. (read only)"""
//|
STATIC mp_obj_t displayio_bitmap_obj_get_width(mp_obj_t self_in) {
    displayio_bitmap_t *self = MP_OBJ_TO_PTR(self_in);

    return MP_OBJ_NEW_SMALL_INT(common_hal_displayio_bitmap_get_width(self));
}

MP_DEFINE_CONST_FUN_OBJ_1(displayio_bitmap_get_width_obj, displayio_bitmap_obj_get_width);

const mp_obj_property_t displayio_bitmap_width_obj = {
    .base.type = &mp_type_property,
    .proxy = {(mp_obj_t)&displayio_bitmap_get_width_obj,
              (mp_obj_t)&mp_const_none_obj,
              (mp_obj_t)&mp_const_none_obj},
};

//|     height: int
//|     """Height of the bitmap. (read only)"""
//|
STATIC mp_obj_t displayio_bitmap_obj_get_height(mp_obj_t self_in) {
    displayio_bitmap_t *self = MP_OBJ_TO_PTR(self_in);

    return MP_OBJ_NEW_SMALL_INT(common_hal_displayio_bitmap_get_height(self));
}

MP_DEFINE_CONST_FUN_OBJ_1(displayio_bitmap_get_height_obj, displayio_bitmap_obj_get_height);

const mp_obj_property_t displayio_bitmap_height_obj = {
    .base.type = &mp_type_property,
    .proxy = {(mp_obj_t)&displayio_bitmap_get_height_obj,
              (mp_obj_t)&mp_const_none_obj,
              (mp_obj_t)&mp_const_none_obj},
};

//|     def __getitem__(self, index: Union[Tuple[int, int], int]) -> int:
//|         """Returns the value at the given index. The index can either be an x,y tuple or an int equal
//|         to ``y * width + x``.
//|
//|         This allows you to::
//|
//|           print(bitmap[0,1])"""
//|         ...
//|
//|     def __setitem__(self, index: Union[Tuple[int, int], int], value: int) -> None:
//|         """Sets the value at the given index. The index can either be an x,y tuple or an int equal
//|         to ``y * width + x``.
//|
//|         This allows you to::
//|
//|           bitmap[0,1] = 3"""
//|         ...
//|
STATIC mp_obj_t bitmap_subscr(mp_obj_t self_in, mp_obj_t index_obj, mp_obj_t value_obj) {
    if (value_obj == mp_const_none) {
        // delete item
        mp_raise_AttributeError(translate("Cannot delete values"));
        return mp_const_none;
    }

    displayio_bitmap_t *self = MP_OBJ_TO_PTR(self_in);

    if (MP_OBJ_IS_TYPE(index_obj, &mp_type_slice)) {
        // TODO(tannewt): Implement subscr after slices support start, stop and step tuples.
        mp_raise_NotImplementedError(translate("Slices not supported"));
        return mp_const_none;
    }

    uint16_t x = 0;
    uint16_t y = 0;
    if (MP_OBJ_IS_SMALL_INT(index_obj)) {
        mp_int_t i = MP_OBJ_SMALL_INT_VALUE(index_obj);
        uint16_t width = common_hal_displayio_bitmap_get_width(self);
        x = i % width;
        y = i / width;
    } else {
        mp_obj_t* items;
        mp_obj_get_array_fixed_n(index_obj, 2, &items);
        x = mp_obj_get_int(items[0]);
        y = mp_obj_get_int(items[1]);
        if (x >= common_hal_displayio_bitmap_get_width(self) || y >= common_hal_displayio_bitmap_get_height(self)) {
            mp_raise_IndexError(translate("pixel coordinates out of bounds"));
        }
    }

    if (value_obj == MP_OBJ_SENTINEL) {
        // load
        return MP_OBJ_NEW_SMALL_INT(common_hal_displayio_bitmap_get_pixel(self, x, y));
    } else {
        mp_uint_t value = (mp_uint_t)mp_obj_get_int(value_obj);
        if ((value >> common_hal_displayio_bitmap_get_bits_per_value(self)) != 0) {
            mp_raise_ValueError(translate("pixel value requires too many bits"));
        }
        common_hal_displayio_bitmap_set_pixel(self, x, y, value);
    }
    return mp_const_none;
}

//|     def blit(self, x: int, y: int, source_bitmap: Bitmap, *, x1: int, y1: int, x2: int, y2: int, skip_index: int) -> None:
//|         """Inserts the source_bitmap region defined by rectangular boundaries
//|                     (x1,y1) and (x2,y2) into the bitmap at the specified (x,y) location.
//|
//|         :param int x: Horizontal pixel location in bitmap where source_bitmap upper-left
//|                       corner will be placed
//|         :param int y: Vertical pixel location in bitmap where source_bitmap upper-left
//|                       corner will be placed
//|         :param bitmap source_bitmap: Source bitmap that contains the graphical region to be copied
//|         :param int x1: Minimum x-value for rectangular bounding box to be copied from the source bitmap
//|         :param int y1: Minimum y-value for rectangular bounding box to be copied from the source bitmap
//|         :param int x2: Maximum x-value (exclusive) for rectangular bounding box to be copied from the source bitmap
//|         :param int y2: Maximum y-value (exclusive) for rectangular bounding box to be copied from the source bitmap
//|         :param int skip_index: bitmap palette index in the source that will not be copied,
//|                                set to None to copy all pixels"""
//|         ...
//|
STATIC mp_obj_t displayio_bitmap_obj_blit(size_t n_args, const mp_obj_t *pos_args, mp_map_t *kw_args){
    enum {ARG_x, ARG_y, ARG_source, ARG_x1, ARG_y1, ARG_x2, ARG_y2, ARG_skip_index};
    static const mp_arg_t allowed_args[] = {
        {MP_QSTR_x, MP_ARG_REQUIRED | MP_ARG_INT},
        {MP_QSTR_y, MP_ARG_REQUIRED | MP_ARG_INT},
        {MP_QSTR_source_bitmap, MP_ARG_REQUIRED | MP_ARG_OBJ},
        {MP_QSTR_x1, MP_ARG_KW_ONLY | MP_ARG_INT, {.u_int = 0} },
        {MP_QSTR_y1, MP_ARG_KW_ONLY | MP_ARG_INT, {.u_int = 0} },
        {MP_QSTR_x2, MP_ARG_KW_ONLY | MP_ARG_OBJ, {.u_obj = mp_const_none} }, // None convert to source->width
        {MP_QSTR_y2, MP_ARG_KW_ONLY | MP_ARG_OBJ, {.u_obj = mp_const_none} }, // None convert to source->height
        {MP_QSTR_skip_index, MP_ARG_OBJ | MP_ARG_KW_ONLY, {.u_obj=mp_const_none} },
    };
    mp_arg_val_t args[MP_ARRAY_SIZE(allowed_args)];
    mp_arg_parse_all(n_args - 1, pos_args + 1, kw_args, MP_ARRAY_SIZE(allowed_args), allowed_args, args);

    displayio_bitmap_t *self = MP_OBJ_TO_PTR(pos_args[0]);

    int16_t x = args[ARG_x].u_int;
    int16_t y = args[ARG_y].u_int;

    displayio_bitmap_t *source = MP_OBJ_TO_PTR(args[ARG_source].u_obj);

    // ensure that the target bitmap (self) has at least as many `bits_per_value` as the source
    if (self->bits_per_value < source->bits_per_value) {
        mp_raise_ValueError(translate("source palette too large"));
    }

    int16_t x1 = args[ARG_x1].u_int;
    int16_t y1 = args[ARG_y1].u_int;
    int16_t x2, y2;
    // if x2 or y2 is None, then set as the maximum size of the source bitmap
    if ( args[ARG_x2].u_obj == mp_const_none ) {
        x2 = source->width;
    } else {
        x2 = mp_obj_get_int(args[ARG_x2].u_obj);
    }
    //int16_t y2;
    if ( args[ARG_y2].u_obj == mp_const_none ) {
        y2 = source->height;
    } else {
        y2 = mp_obj_get_int(args[ARG_y2].u_obj);
    }

    // Check x,y are within self (target) bitmap boundary
    if ( (x < 0) || (y < 0) || (x > self->width) || (y > self->height) ) {
            mp_raise_ValueError(translate("out of range of target"));
    }
    // Check x1,y1,x2,y2 are within source bitmap boundary
    if ( (x1 < 0) || (x1 > source->width)  ||
        (y1 < 0) || (y1 > source->height) ||
        (x2 < 0) || (x2 > source->width)  ||
        (y2 < 0) || (y2 > source->height) ) {
            mp_raise_ValueError(translate("out of range of source"));
    }

    // Ensure x1 < x2 and y1 < y2
    if (x1 > x2) {
        int16_t temp=x2;
        x2=x1;
        x1=temp;
    }
    if (y1 > y2) {
        int16_t temp=y2;
        y2=y1;
        y1=temp;
    }

    uint32_t skip_index;
    bool skip_index_none; // flag whether skip_value was None

    if (args[ARG_skip_index].u_obj == mp_const_none ) {
        skip_index = 0;
        skip_index_none = true;
    } else {
        skip_index = mp_obj_get_int(args[ARG_skip_index].u_obj);
        skip_index_none = false;
    }

    common_hal_displayio_bitmap_blit(self, x, y, source, x1, y1, x2, y2, skip_index, skip_index_none);

    return mp_const_none;
}
MP_DEFINE_CONST_FUN_OBJ_KW(displayio_bitmap_blit_obj, 4, displayio_bitmap_obj_blit);
// `displayio_bitmap_obj_blit` requires at least 4 arguments

/*
STATIC int16_t validate_point(mp_obj_t point, int16_t default_value) {
    // Checks if point is None and returns default_value, otherwise decodes integer value
    if ( point == mp_const_none ) {
        return default_value;
    }
    return mp_obj_get_int(point);
}

STATIC void extract_tuple(mp_obj_t xy_tuple, int16_t *x, int16_t *y, int16_t x_default, int16_t y_default) {
    // Helper function for fancyblit
    // Extract x,y values from a tuple or default if None
    if ( xy_tuple == mp_const_none ) {
        *x = x_default;
        *y = y_default;
    } else if ( !MP_OBJ_IS_OBJ(xy_tuple) ) {
        mp_raise_ValueError(translate("clip point must be (x,y) tuple"));
    } else {
        mp_obj_t* items;
        mp_obj_get_array_fixed_n(xy_tuple, 2, &items);
        *x = mp_obj_get_int(items[0]);
        *y = mp_obj_get_int(items[1]);
    }
}

STATIC void validate_clip_region(displayio_bitmap_t *bitmap, mp_obj_t clip0_tuple, int16_t *clip0_x, int16_t *clip0_y,
                                                   mp_obj_t clip1_tuple, int16_t *clip1_x, int16_t *clip1_y) {
    // Helper function for fancyblit
    // 1. Extract the clip x,y points from the two clip tuples
    // 2. Rearrange values such that clip0_ < clip1_
    // 3. Constrain the clip points to within the bitmap

    extract_tuple(clip0_tuple, clip0_x, clip0_y, 0, 0);
    extract_tuple(clip1_tuple, clip1_x, clip1_y, bitmap->width, bitmap->height);

    // Ensure the value for clip0 is less than clip1 (for both x and y)
    if ( *clip0_x > *clip1_x ) {
        int16_t temp_value = *clip0_x; // swap values
        *clip0_x = *clip1_x;
        *clip1_x = temp_value;
    }
    if ( *clip0_y > *clip1_y ) {
        int16_t temp_value = *clip0_y; // swap values
        *clip0_y = *clip1_y;
        *clip1_y = temp_value;
    }

    // Constrain the clip window to within the bitmap boundaries
    if (*clip0_x < 0) {
        *clip0_x = 0;
    }
    if (*clip0_y < 0) {
        *clip0_y = 0;
    }
    if (*clip0_x > bitmap->width) {
        *clip0_x = bitmap->width;
    }
    if (*clip0_y > bitmap->height) {
        *clip0_y = bitmap->height;
    }
    if (*clip1_x < 0) {
        *clip1_x = 0;
    }
    if (*clip1_y < 0) {
        *clip1_y = 0;
    }
    if (*clip1_x > bitmap->width) {
        *clip1_x = bitmap->width;
    }
    if (*clip1_y > bitmap->height) {
        *clip1_y = bitmap->height;
    }

}

/////////  A fancy blit "rotozoom" function with rotation, scaling and clipping (both source and destination)
//|     def rotozoom(self, ox: int, oy: int,
//|                       dest_clip0: Tuple[int, int],
//|                       dest_clip1: Tuple[int, int],
//|                       source_bitmap: Bitmap,
//|                       px: int, py: int,
//|                       source_clip0: Tuple[int, int],
//|                       source_clip1: Tuple[int, int],
//|                       angle: float,
//|                       scale: float,
//|                       skip_index: int) -> None:
//|         """Inserts the source_bitmap region into the bitmap with rotation (angle), scale
//|                       and clipping (both on source and destination bitmaps)
//|
//|         :param int ox: Horizontal pixel location in destination bitmap where source bitmap
//|                       point (px,py) is placed
//|         :param int oy: Vertical pixel location in destination bitmap where source bitmap
//|                       point (px,py) is placed
//|         :param dest_clip0: First corner of rectangular destination clipping
//|                       region that constrains region of writing into destination bitmap
//|         :type dest_clip0: Tuple[int,int]
//|         :param dest_clip1: second corner of rectangular destination clipping
//|                       region that constrains region of writing into destination bitmap
//|         :type dest_clip1: Tuple[int,int]
//|         :param bitmap source_bitmap: Source bitmap that contains the graphical region to be copied
//|         :param int px: Horizontal pixel location in source bitmap that is placed into the
//|                       destination bitmap at (ox,oy)
//|         :param int py: Vertical pixel location in source bitmap that is placed into the
//|                       destination bitmap at (ox,oy)
//|         :param source_clip0: First corner of rectangular destination clipping
//|                       region that constrains region of writing into destination bitmap
//|         :type source_clip0: Tuple[int,int]
//|         :param source_clip1: second corner of rectangular destination clipping
//|                       region that constrains region of writing into destination bitmap
//|         :type source_clip1: Tuple[int,int]
//|         :param float angle: angle of rotation, in radians (positive is clockwise direction)
//|         :param float scale: scaling factor
//|         :param int skip_index: bitmap palette index in the source that will not be copied,
//|                                set to None to copy all pixels"""
//|         ...
//|
STATIC mp_obj_t displayio_bitmap_obj_rotozoom(size_t n_args, const mp_obj_t *pos_args, mp_map_t *kw_args){
    enum {ARG_ox, ARG_oy, ARG_dest_clip0, ARG_dest_clip1,
        ARG_source_bitmap, ARG_px, ARG_py,
        ARG_source_clip0, ARG_source_clip1,
        ARG_angle, ARG_scale, ARG_skip_index};

    static const mp_arg_t allowed_args[] = {
        {MP_QSTR_ox, MP_ARG_KW_ONLY | MP_ARG_OBJ, {.u_obj = mp_const_none} }, // None convert to destination->width  / 2
        {MP_QSTR_oy, MP_ARG_KW_ONLY | MP_ARG_OBJ, {.u_obj = mp_const_none} }, // None convert to destination->height / 2
        {MP_QSTR_dest_clip0, MP_ARG_KW_ONLY | MP_ARG_OBJ, {.u_obj = mp_const_none} },
        {MP_QSTR_dest_clip1, MP_ARG_KW_ONLY | MP_ARG_OBJ, {.u_obj = mp_const_none} },

        {MP_QSTR_source_bitmap, MP_ARG_REQUIRED | MP_ARG_OBJ},
        {MP_QSTR_px, MP_ARG_KW_ONLY | MP_ARG_OBJ, {.u_obj = mp_const_none} }, // None convert to source->width  / 2
        {MP_QSTR_py, MP_ARG_KW_ONLY | MP_ARG_OBJ, {.u_obj = mp_const_none} }, // None convert to source->height / 2
        {MP_QSTR_source_clip0, MP_ARG_KW_ONLY | MP_ARG_OBJ, {.u_obj = mp_const_none} },
        {MP_QSTR_source_clip1, MP_ARG_KW_ONLY | MP_ARG_OBJ, {.u_obj = mp_const_none} },

        {MP_QSTR_angle, MP_ARG_KW_ONLY | MP_ARG_OBJ, {.u_obj = mp_const_none} }, // None convert to 0.0
        {MP_QSTR_scale, MP_ARG_KW_ONLY | MP_ARG_OBJ, {.u_obj = mp_const_none} }, // None convert to 1.0
        {MP_QSTR_skip_index, MP_ARG_OBJ | MP_ARG_KW_ONLY, {.u_obj=mp_const_none} },
    };

    mp_arg_val_t args[MP_ARRAY_SIZE(allowed_args)];
    mp_arg_parse_all(n_args - 1, pos_args + 1, kw_args, MP_ARRAY_SIZE(allowed_args), allowed_args, args);

    displayio_bitmap_t *self = MP_OBJ_TO_PTR(pos_args[0]); // self is the destination bitmap

    displayio_bitmap_t *source = MP_OBJ_TO_PTR(args[ARG_source_bitmap].u_obj); // the source bitmap

    // ensure that the target bitmap (self) has at least as many `bits_per_value` as the source
    if (self->bits_per_value < source->bits_per_value) {
        mp_raise_ValueError(translate("source palette too large"));
    }

    // Confirm the destination location target (ox,oy); if None, default to bitmap midpoint
    int16_t ox, oy;
    ox = validate_point(args[ARG_ox].u_obj, self->width  / 2);
    oy = validate_point(args[ARG_oy].u_obj, self->height / 2);

    // Confirm the source location target (px,py); if None, default to bitmap midpoint
    int16_t px, py;
    px = validate_point(args[ARG_px].u_obj, source->width  / 2);
    py = validate_point(args[ARG_py].u_obj, source->height / 2);

    // Validate the clipping regions for the destination bitmap
    int16_t dest_clip0_x, dest_clip0_y, dest_clip1_x, dest_clip1_y;

    validate_clip_region(self, args[ARG_dest_clip0].u_obj, &dest_clip0_x, &dest_clip0_y,
                               args[ARG_dest_clip1].u_obj, &dest_clip1_x, &dest_clip1_y);

    // Validate the clipping regions for the source bitmap
    int16_t source_clip0_x, source_clip0_y, source_clip1_x, source_clip1_y;

    validate_clip_region(source, args[ARG_source_clip0].u_obj, &source_clip0_x, &source_clip0_y,
                                 args[ARG_source_clip1].u_obj, &source_clip1_x, &source_clip1_y);

    // Confirm the angle value
    float angle=0.0;
    if ( args[ARG_angle].u_obj != mp_const_none ) {
        angle = mp_obj_get_float(args[ARG_angle].u_obj);
    }

    // Confirm the scale value
    float scale=1.0;
    if ( args[ARG_scale].u_obj != mp_const_none ) {
        scale = mp_obj_get_float(args[ARG_scale].u_obj);
    }
    if (scale < 0) { // ensure scale >= 0
        scale = 1.0;
    }

    uint32_t skip_index;
    bool skip_index_none; // Flag whether input skip_value was None
    if (args[ARG_skip_index].u_obj == mp_const_none ) {
        skip_index = 0;
        skip_index_none = true;
    } else {
        skip_index = mp_obj_get_int(args[ARG_skip_index].u_obj);
        skip_index_none = false;
    }

    common_hal_displayio_bitmap_rotozoom(self, ox, oy,
                                            dest_clip0_x, dest_clip0_y,
                                            dest_clip1_x, dest_clip1_y,
                                            source, px, py,
                                            source_clip0_x, source_clip0_y,
                                            source_clip1_x, source_clip1_y,
                                            angle,
                                            scale,
                                            skip_index, skip_index_none);

    return mp_const_none;
}

MP_DEFINE_CONST_FUN_OBJ_KW(displayio_bitmap_rotozoom_obj, 1, displayio_bitmap_obj_rotozoom);
// requires at least 2 arguments (destination and source bitmaps)
*/


//|     def fill(self, value: int) -> None:
//|         """Fills the bitmap with the supplied palette index value."""
//|         ...
//|
STATIC mp_obj_t displayio_bitmap_obj_fill(mp_obj_t self_in, mp_obj_t value_obj) {
    displayio_bitmap_t *self = MP_OBJ_TO_PTR(self_in);

    mp_uint_t value = (mp_uint_t)mp_obj_get_int(value_obj);
    if ((value >> common_hal_displayio_bitmap_get_bits_per_value(self)) != 0) {
            mp_raise_ValueError(translate("pixel value requires too many bits"));
    }
    common_hal_displayio_bitmap_fill(self, value);

    return mp_const_none;
}
MP_DEFINE_CONST_FUN_OBJ_2(displayio_bitmap_fill_obj, displayio_bitmap_obj_fill);

STATIC const mp_rom_map_elem_t displayio_bitmap_locals_dict_table[] = {
    { MP_ROM_QSTR(MP_QSTR_height), MP_ROM_PTR(&displayio_bitmap_height_obj) },
    { MP_ROM_QSTR(MP_QSTR_width), MP_ROM_PTR(&displayio_bitmap_width_obj) },
    { MP_ROM_QSTR(MP_QSTR_blit), MP_ROM_PTR(&displayio_bitmap_blit_obj) },
//    { MP_ROM_QSTR(MP_QSTR_rotozoom), MP_ROM_PTR(&displayio_bitmap_rotozoom_obj) },
    { MP_ROM_QSTR(MP_QSTR_fill), MP_ROM_PTR(&displayio_bitmap_fill_obj) },

};
STATIC MP_DEFINE_CONST_DICT(displayio_bitmap_locals_dict, displayio_bitmap_locals_dict_table);

const mp_obj_type_t displayio_bitmap_type = {
    { &mp_type_type },
    .name = MP_QSTR_Bitmap,
    .make_new = displayio_bitmap_make_new,
    .subscr = bitmap_subscr,
    .locals_dict = (mp_obj_dict_t*)&displayio_bitmap_locals_dict,
};
