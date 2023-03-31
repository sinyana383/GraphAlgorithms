#ifndef NAVIGATOR_SRC_STACK_HPP_
#define NAVIGATOR_SRC_STACK_HPP_

#include <cstdlib>

template <typename T>  // НЕ адаптированно под классы
struct s_list {
  T content;
  struct s_list *next = nullptr;

  s_list() = default;
  s_list(T c, struct s_list *n) : content(c), next(n) {}
};
template <typename T>
using t_list = struct s_list<T>;  // alias template - name that refers to a
                                  // previously defined type

template <typename T>
class Stack {
 private:
  t_list<T> *_top = nullptr;
  std::size_t size = 0;

 public:
  void init() { _top = new t_list<T>; }
  void push(T value);
  T pop();
  T peek() {
    return _top->next != nullptr ? _top->next->content : _top->content;
  }
  ~Stack();

  size_t GetSize() const;
};

template <typename T>
void Stack<T>::push(T value) {
  auto *newbie = new t_list<T>(value, _top->next);
  _top->next = newbie;
  ++size;
}

template <typename T>
T Stack<T>::pop() {
  T value = this->peek();
  auto todel = _top->next;

  _top->next = todel != nullptr ? todel->next : nullptr;
  delete todel;

  if (size == 0) throw std::out_of_range("pop from empty stack");
  --size;
  return value;
}
template <typename T>
size_t Stack<T>::GetSize() const {
  return size;
}
template <typename T>
Stack<T>::~Stack() {
  while (_top->next != nullptr) this->pop();
  delete[] _top;
}

#endif  // NAVIGATOR_SRC_STACK_HPP_
