#pragma once

#include <amal/internal/vec2.hpp>
#include <amal/internal/vec3.hpp>
#include <amal/internal/vec4.hpp>
#include <amal/internal/vector/operators.hpp>

#ifdef AMAL_INTEGRATE_ACUL
    #include <acul/string/utils.hpp>

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
} // namespace acul
#endif