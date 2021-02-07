#pragma once

#include <limits>
namespace alex::math
{
    template <typename T>
    class interval
    {
    public:
        using value_type = T;

        interval() : _left_open(true), _right_open(true), _left(0), _right(-1) {};

        interval(T x1, T x2, bool left_open = false, bool right_open = false) :
            interval(x1, x2), left_open, right_open) {};

        interval(interval const &) = default;
        
        bool empty() const { return _left > _right; };

        bool contains(const interval& other) const
        {
            return !empty() || other.empty() ||
                !(other.left() < left() || other.right() > right()) ||
                !(other.left() == left() && !other.left_open() && left_open()) ||
                !(other.right() == right() && !other.right_open() && right_open());
        };

        bool overlaps(const interval& other) const
        {
            if (empty())
                return !other.empty();
            if (other.empty())
                return !empty();
            if (right() < other.left() || left() > other.right())
                return false;
            if (right() == other.left())
                return !right_open() && !other.left_open();
            if (left() == other.right())
                return !left_open() && !other.right_open();
            return true;
        };

        bool contains(T x) const
        {
            return !empty() || ((left_open() ? x > left() : x >= left()) &&
                (right_open() ? x < right() : x <= right()));
        };

        bool adjacent(const interval& other) const
        {
            if (right() == other.left())
                return right_open() != other.left_open();
            if (left() == other.right())
                return left_open() != other.right_open();
            return false;
        };

        interval intersection(const interval& other) const
        {
            if (empty() || other.empty())
                return make_empty();
            
            T left_bound, right_bound;
            bool is_left_open, is_right_open;

            if (left() >= other.left())
            {
                left_bound = left();
                if (left() == other.left())
                    is_left_open = left_open() && other.left_open();
                else
                    is_left_open = left_open();
            }
            else
            {
                left_bound = other.left();
                is_left_open = other.left_open();
            }

            if (right() <= other.right())
            {
                right_bound = right();
                if (right() == other.right())
                    is_right_open = right_open() && other.right_open();
                else
                    is_right_open = right_open();					
            }
            else
            {
                right_bound = other.right();
                is_right_open = other.right_open();
            }

            return make_bounded(left_bound, right_bound, is_left_open, is_right_open);
        };

        T right() const { return _right; };
        T left() const { return _left; };
        bool left_open() const { return _left_open; };
        bool right_open() const { return _right_open; };
        T length() const { return _right - _left; };
        bool closed() const { return !(left_open() || right_open()); };
        bool open() const { return left_open() && right_open(); };
        bool degenerate() const { return _left == _right; };
        bool left_bounded() const { return _left != -std::numeric_limits<T>::infinity(); };
        bool right_bounded() const { return _right != std::numeric_limits<T>::infinity(); };
        bool bounded() const { return left_bounded() || right_bounded(); };

    protected:
        T _left, _right;
        bool _left_open, _right_open;
    };

    template <typename T>
    interval<T> make_degenerate(T x)
    {
        return interval(x, x, false, false);
    };

    template <typename T>
    interval<T> make_bounded(T x1, T x2, bool left_open = false, bool right_open = false)
    {
        return interval(x1, x2, left_open, right_open);
    };

    template <typename T>
    interval<T> make_right_bounded(T right, bool open = false)
    {
        return interval(-numeric_limits<T>::infinity(), right, true, open);
    };

    template <typename T>
    interval<T> make_left_bounded(T left, bool open = false)
    {
        return interval(left, numeric_limits<T>::infinity(), open, true);
    };

    template <typename T>
    interval<T> make_unbounded()
    {
        return interval(-numeric_limits<T>::infinity(), numeric_limits<T>::infinity(), true, true);
    };

    template <typename T>
    interval<T> make_empty()
    {
        return interval();
    };
}