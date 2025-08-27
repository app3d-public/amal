#pragma once

#include "../../vector.hpp"

namespace acul
{
    template <>
    bin_stream &write(const amal::vec2 &vec)
    {
        return write(vec.x).write(vec.y);
    }

    template <>
    bin_stream &write(const amal::vec3 &vec)
    {
        return write(vec.x).write(vec.y).write(vec.z);
    }

    template <>
    bin_stream &write(const amal::vec4 &vec)
    {
        return write(vec.x).write(vec.y).write(vec.z).write(vec.w);
    }

    template <>
    bin_stream &write(const amal::ivec2 &vec)
    {
        return write(vec.x).write(vec.y);
    }

    template <>
    bin_stream &write(const amal::ivec3 &vec)
    {
        return write(vec.x).write(vec.y).write(vec.z);
    }

    template <>
    bin_stream &write(const amal::ivec4 &vec)
    {
        return write(vec.x).write(vec.y).write(vec.z).write(vec.w);
    }

    template <>
    bin_stream &write(const amal::dvec2 &vec)
    {
        return write(vec.x).write(vec.y);
    }

    template <>
    bin_stream &write(const amal::dvec3 &vec)
    {
        return write(vec.x).write(vec.y).write(vec.z);
    }

    template <>
    bin_stream &write(const amal::dvec4 &vec)
    {
        return write(vec.x).write(vec.y).write(vec.z).write(vec.w);
    }

    template <>
    bin_stream &write(const amal::mat2 &mat)
    {
        return write(mat[0][0]).write(mat[0][1]).write(mat[1][0]).write(mat[1][1]);
    }

    template <>
    bin_stream &write(const amal::mat3 &mat)
    {
        return write(mat[0][0])
            .write(mat[0][1])
            .write(mat[0][2])
            .write(mat[1][0])
            .write(mat[1][1])
            .write(mat[1][2])
            .write(mat[2][0])
            .write(mat[2][1])
            .write(mat[2][2]);
    }

    template <>
    bin_stream &write(const amal::mat4 &mat)
    {
        return write(mat[0][0])
            .write(mat[0][1])
            .write(mat[0][2])
            .write(mat[0][3])
            .write(mat[1][0])
            .write(mat[1][1])
            .write(mat[1][2])
            .write(mat[1][3])
            .write(mat[2][0])
            .write(mat[2][1])
            .write(mat[2][2])
            .write(mat[2][3])
            .write(mat[3][0])
            .write(mat[3][1])
            .write(mat[3][2])
            .write(mat[3][3]);
    }

    template <>
    bin_stream &read(amal::vec2 &vec)
    {
        return read(vec.x).read(vec.y);
    }

    template <>
    bin_stream &read(amal::vec3 &vec)
    {
        return read(vec.x).read(vec.y).read(vec.z);
    }

    template <>
    bin_stream &read(amal::vec4 &vec)
    {
        return read(vec.x).read(vec.y).read(vec.z).read(vec.w);
    }

    template <>
    bin_stream &read(amal::ivec2 &vec)
    {
        return read(vec.x).read(vec.y);
    }

    template <>
    bin_stream &read(amal::ivec3 &vec)
    {
        return read(vec.x).read(vec.y).read(vec.z);
    }

    template <>
    bin_stream &read(amal::ivec4 &vec)
    {
        return read(vec.x).read(vec.y).read(vec.z).read(vec.w);
    }

    template <>
    bin_stream &read(amal::dvec2 &vec)
    {
        return read(vec.x).read(vec.y);
    }

    template <>
    bin_stream &read(amal::dvec3 &vec)
    {
        return read(vec.x).read(vec.y).read(vec.z);
    }

    template <>
    bin_stream &read(amal::dvec4 &vec)
    {
        return read(vec.x).read(vec.y).read(vec.z).read(vec.w);
    }

    bin_stream &read(amal::mat2 &mat) { return read(mat[0][0]).read(mat[0][1]).read(mat[1][0]).read(mat[1][1]); }

    bin_stream &read(amal::mat3 &mat)
    {
        return read(mat[0][0])
            .read(mat[0][1])
            .read(mat[0][2])
            .read(mat[1][0])
            .read(mat[1][1])
            .read(mat[1][2])
            .read(mat[2][0])
            .read(mat[2][1])
            .read(mat[2][2]);
    }

    bin_stream &read(amal::mat4 &mat)
    {
        return read(mat[0][0])
            .read(mat[0][1])
            .read(mat[0][2])
            .read(mat[0][3])
            .read(mat[1][0])
            .read(mat[1][1])
            .read(mat[1][2])
            .read(mat[1][3])
            .read(mat[2][0])
            .read(mat[2][1])
            .read(mat[2][2])
            .read(mat[2][3])
            .read(mat[3][0])
            .read(mat[3][1])
            .read(mat[3][2])
            .read(mat[3][3]);
    }
} // namespace acul