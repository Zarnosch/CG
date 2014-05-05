#include "PerspectiveCamera.hpp"

namespace rt {

Ray PerspectiveCamera::ray(size_t x, size_t y) const
{
  return Ray(this->position(),
           ((this->topLeft() + this->right()*double(x) - this->down()*double(y)) - this->position()));
}

} //namespace rt
