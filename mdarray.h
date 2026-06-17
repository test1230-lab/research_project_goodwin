#include <vector>
#include <cstdlib>
#pragma once

template<typename T>
class array3d
{
public:
    array3d() = default;
    array3d(std::size_t N, std::size_t M, std::size_t P) 
        : vec(N*M*P), N(N), M(M), P(P){}
    array3d(std::size_t N, std::size_t M, std::size_t P, T val) 
        : vec(N*M*P, val), N(N), M(M), P(P){}

    T& operator[](std::size_t i, std::size_t j, std::size_t k) noexcept
    {
        return vec[(i*M + j)*P + k];
    }

    const T& operator[](std::size_t i, std::size_t j, std::size_t k) const noexcept
    {
        return vec[(i*M + j)*P + k];
    }

    void resize(std::size_t N, std::size_t M, std::size_t P)
    {
        this->N = N;
        this->M = M;
        this->P = P;
        vec.resize(N*M*P);
    }

    T* data() noexcept { return vec.data(); }
    const T* data() const noexcept { return vec.data(); }

    std::size_t dim0() const noexcept { return N; }
    std::size_t dim1() const noexcept { return M; }
    std::size_t dim2() const noexcept { return P; }

    std::size_t size() const noexcept { return vec.size(); }

    auto begin() noexcept { return vec.begin(); }
    auto end() noexcept { return vec.end(); }

    auto begin() const noexcept { return vec.begin(); }
    auto end() const noexcept { return vec.end(); }

    auto cbegin() const noexcept { return vec.cbegin(); }
    auto cend() const noexcept { return vec.cend(); }

private:
    std::size_t N, M, P;
    std::vector<T> vec;
};

template<typename T>
class array2d
{
public:
    array2d() = default;
    array2d(std::size_t N, std::size_t M) 
        : vec(N*M), N(N), M(M){}
    array2d(std::size_t N, std::size_t M, T val) 
        : vec(N*M, val), N(N), M(M){}

    T& operator[](std::size_t i, std::size_t j) noexcept
    {
        return vec[i*M + j];
    }

    const T& operator[](std::size_t i, std::size_t j) const noexcept
    {
        return vec[i*M + j];
    }

    void resize(std::size_t N, std::size_t M)
    {
        this->N = N;
        this->M = M;
        vec.resize(N*M);
    }

    T* data() noexcept { return vec.data(); }
    const T* data() const noexcept { return vec.data(); }

    std::size_t dim0() const noexcept { return N; }
    std::size_t dim1() const noexcept { return M; }

    std::size_t size() const noexcept { return vec.size(); }

    auto begin() noexcept { return vec.begin(); }
    auto end() noexcept { return vec.end(); }

    auto begin() const noexcept { return vec.begin(); }
    auto end() const noexcept { return vec.end(); }

    auto cbegin() const noexcept { return vec.cbegin(); }
    auto cend() const noexcept { return vec.cend(); }

private:
    std::size_t N, M;
    std::vector<T> vec;
};