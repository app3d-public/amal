#pragma once

#include <acul/string/utils.hpp>
#include "../../matrix.hpp"
#include "../../vector.hpp"

namespace acul
{
    /**
     * @brief Converts the 2-dimensional vector to the C-style string
     * @param vec Source vector
     * @param buffer Destination buffer
     * @param buffer_size Buffer size
     * @return The number of characters written
     **/
    inline int to_string(const amal::vec2 &vec, char *buffer, size_t buffer_size, size_t offset)
    {
        int written = to_string(vec.x, buffer + offset, buffer_size - offset, 5);
        if (written == 0) return 0;
        offset += written;
        if (offset < buffer_size - 1)
            buffer[offset++] = ' ';
        else
            return 0;

        written = to_string(vec.y, buffer + offset, buffer_size - offset, 5);
        if (written == 0) return 0;
        offset += written;

        return offset;
    }

    /**
     * @brief Converts the 3-dimensional vector to the C-style string
     * @param vec Source vector
     * @param buffer Destination buffer
     * @param buffer_size Buffer size
     * @return The number of characters written
     **/
    inline int to_string(const amal::vec3 &vec, char *buffer, size_t buffer_size, size_t offset)
    {
        int written = to_string(vec.x, buffer + offset, buffer_size - offset, 5);
        if (written == 0) return 0;
        offset += written;
        if (offset < buffer_size - 1)
            buffer[offset++] = ' ';
        else
            return 0;

        written = to_string(vec.y, buffer + offset, buffer_size - offset, 5);
        if (written == 0) return 0;
        offset += written;
        if (offset < buffer_size - 1)
            buffer[offset++] = ' ';
        else
            return 0;

        written = to_string(vec.z, buffer + offset, buffer_size - offset, 5);
        if (written == 0) return 0;
        offset += written;

        return offset;
    }

    /**
     * @brief Deserialize the C-style string to the 2-dimensional vector.
     * All values to deserialize must be present in the string.
     * @param str Source string
     * @param vec Destination vector
     * @return True if successful. Otherwise false
     **/
    inline bool stov2(const char *&str, amal::vec2 &vec) { return stof(str, vec.x) && stof(str, vec.y); }

    /**
     * @brief Deserialize the C-style string to the 2-dimensional vector.
     * This function deserializes line, attempting to fill as many components as possible based on the input.
     * @param str Source string
     * @param vec Destination vector
     * @return True if successful. Otherwise false
     **/
    inline void stov2_opt(const char *&str, amal::vec2 &vec)
    {
        if (!stof(str, vec.x)) return;
        if (!stof(str, vec.y)) return;
    }

    /**
     * @brief Deserialize the C-style string to the 3-dimensional vector.
     * All values to deserialize must be present in the string.
     * @param str Source string
     * @param vec Destination vector
     * @return True if successful. Otherwise false
     **/
    inline bool stov3(const char *&str, amal::vec3 &vec)
    {
        return stof(str, vec.x) && stof(str, vec.y) && stof(str, vec.z);
    }

    /**
     * @brief Deserialize the C-style string to the 3-dimensional vector.
     * This function deserializes line, attempting to fill as many components as possible based on the input.
     * @param str Source string
     * @param vec Destination vector
     * @return True if successful. Otherwise false
     **/
    inline void stov3_opt(const char *&str, amal::vec3 &vec)
    {
        if (!stof(str, vec.x)) return;
        if (!stof(str, vec.y)) return;
        if (!stof(str, vec.z)) return;
    }

    inline string to_string(const amal::bvec2 &v) { return format("bvec2(%d, %d)", v.x, v.y); }

    inline string to_string(const amal::bvec3 &v) { return format("bvec3(%d, %d, %d)", v.x, v.y, v.z); }

    inline string to_string(const amal::bvec4 &v) { return format("bvec4(%d, %d, %d, %d)", v.x, v.y, v.z, v.w); }

    inline string to_string(const amal::ivec2 &v) { return format("ivec2(%d, %d)", v.x, v.y); }

    inline string to_string(const amal::ivec3 &v) { return format("ivec3(%d, %d, %d)", v.x, v.y, v.z); }

    inline string to_string(const amal::ivec4 &v) { return format("ivec4(%d, %d, %d, %d)", v.x, v.y, v.z, v.w); }

    inline string to_string(const amal::vec2 &v) { return format("vec2(%f, %f)", v.x, v.y); }

    inline string to_string(const amal::vec3 &v) { return format("vec3(%f, %f, %f)", v.x, v.y, v.z); }

    inline string to_string(const amal::vec4 &v) { return format("vec4(%f, %f, %f, %f)", v.x, v.y, v.z, v.w); }

    inline string to_string(const amal::dvec2 &v) { return format("dvec2(%lf, %lf)", v.x, v.y); }

    inline string to_string(const amal::dvec3 &v) { return format("dvec3(%lf, %lf, %lf)", v.x, v.y, v.z); }

    inline string to_string(const amal::dvec4 &v) { return format("dvec4(%lf, %lf, %lf, %lf)", v.x, v.y, v.z, v.w); }

    inline string to_string(const amal::mat2 &m)
    {
        return format("mat2([[%f, %f]; [%f, %f]])", m[0].x, m[1].x, m[0].y, m[1].y);
    }

    inline string to_string(const amal::mat3 &m)
    {
        return format("mat3([[%f, %f, %f]; [%f, %f, %f]; [%f, %f, %f]])", m[0].x, m[1].x, m[2].x, m[0].y, m[1].y,
                      m[2].y, m[0].z, m[1].z, m[2].z);
    }

    inline string to_string(const amal::mat4 &m)
    {
        return format("mat4([[%f, %f, %f, %f]; [%f, %f, %f, %f]; [%f, %f, %f, %f]; [%f, %f, %f, %f]])", m[0].x, m[1].x,
                      m[2].x, m[3].x, m[0].y, m[1].y, m[2].y, m[3].y, m[0].z, m[1].z, m[2].z, m[3].z, m[0].w, m[1].w,
                      m[2].w, m[3].w);
    }

    inline string to_string(const amal::mat2x3 &m)
    {
        return format("mat2x3([[%f, %f]; [%f, %f]; [%f, %f]])", m[0].x, m[1].x, m[0].y, m[1].y, m[0].z, m[1].z);
    }

    inline string to_string(const amal::mat3x2 &m)
    {
        return format("mat3x2([[%f, %f, %f]; [%f, %f, %f]])", m[0].x, m[1].x, m[2].x, m[0].y, m[1].y, m[2].y);
    }

    inline string to_string(const amal::mat2x4 &m)
    {
        return format("mat2x4([[%f, %f]; [%f, %f]; [%f, %f]; [%f, %f]])", m[0].x, m[1].x, m[0].y, m[1].y, m[0].z,
                      m[1].z, m[0].w, m[1].w);
    }

    inline string to_string(const amal::mat4x2 &m)
    {
        return format("mat4x2([[%f, %f, %f, %f]; [%f, %f, %f, %f]])", m[0].x, m[1].x, m[2].x, m[3].x, m[0].y, m[1].y,
                      m[2].y, m[3].y);
    }

    inline string to_string(const amal::mat3x4 &m)
    {
        return format("mat3x4([[%f, %f, %f]; [%f, %f, %f]; [%f, %f, %f]; [%f, %f, %f]])", m[0].x, m[1].x, m[2].x,
                      m[0].y, m[1].y, m[2].y, m[0].z, m[1].z, m[2].z, m[0].w, m[1].w, m[2].w);
    }

    inline string to_string(const amal::mat4x3 &m)
    {
        return format("mat4x3([[%f, %f, %f, %f]; [%f, %f, %f, %f]; [%f, %f, %f, %f]])", m[0].x, m[1].x, m[2].x, m[3].x,
                      m[0].y, m[1].y, m[2].y, m[3].y, m[0].z, m[1].z, m[2].z, m[3].z);
    }

    inline string to_string(const amal::dmat2 &m)
    {
        return format("dmat2([[%lf, %lf]; [%lf, %lf]])", m[0].x, m[1].x, m[0].y, m[1].y);
    }
    inline string to_string(const amal::dmat3 &m)
    {
        return format("dmat3([[%lf, %lf, %lf]; [%lf, %lf, %lf]; [%lf, %lf, %lf]])", m[0].x, m[1].x, m[2].x, m[0].y,
                      m[1].y, m[2].y, m[0].z, m[1].z, m[2].z);
    }
    inline string to_string(const amal::dmat4 &m)
    {
        return format("dmat4([[%lf, %lf, %lf, %lf]; [%lf, %lf, %lf, %lf]; [%lf, %lf, %lf, %lf]; [%lf, %lf, %lf, %lf]])",
                      m[0].x, m[1].x, m[2].x, m[3].x, m[0].y, m[1].y, m[2].y, m[3].y, m[0].z, m[1].z, m[2].z, m[3].z,
                      m[0].w, m[1].w, m[2].w, m[3].w);
    }

    inline string to_string(const amal::dmat2x3 &m)
    {
        return format("dmat2x3([[%lf, %lf]; [%lf, %lf]; [%lf, %lf]])", m[0].x, m[1].x, m[0].y, m[1].y, m[0].z, m[1].z);
    }
    inline string to_string(const amal::dmat3x2 &m)
    {
        return format("dmat3x2([[%lf, %lf, %lf]; [%lf, %lf, %lf]])", m[0].x, m[1].x, m[2].x, m[0].y, m[1].y, m[2].y);
    }
    inline string to_string(const amal::dmat2x4 &m)
    {
        return format("dmat2x4([[%lf, %lf]; [%lf, %lf]; [%lf, %lf]; [%lf, %lf]])", m[0].x, m[1].x, m[0].y, m[1].y,
                      m[0].z, m[1].z, m[0].w, m[1].w);
    }
    inline string to_string(const amal::dmat4x2 &m)
    {
        return format("dmat4x2([[%lf, %lf, %lf, %lf]; [%lf, %lf, %lf, %lf]])", m[0].x, m[1].x, m[2].x, m[3].x, m[0].y,
                      m[1].y, m[2].y, m[3].y);
    }
    inline string to_string(const amal::dmat3x4 &m)
    {
        return format("dmat3x4([[%lf, %lf, %lf]; [%lf, %lf, %lf]; [%lf, %lf, %lf]; [%lf, %lf, %lf]])", m[0].x, m[1].x,
                      m[2].x, m[0].y, m[1].y, m[2].y, m[0].z, m[1].z, m[2].z, m[0].w, m[1].w, m[2].w);
    }
    inline string to_string(const amal::dmat4x3 &m)
    {
        return format("dmat4x3([[%lf, %lf, %lf, %lf]; [%lf, %lf, %lf, %lf]; [%lf, %lf, %lf, %lf]])", m[0].x, m[1].x,
                      m[2].x, m[3].x, m[0].y, m[1].y, m[2].y, m[3].y, m[0].z, m[1].z, m[2].z, m[3].z);
    }
} // namespace acul