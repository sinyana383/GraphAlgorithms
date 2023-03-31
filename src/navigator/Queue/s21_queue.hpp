#ifndef NAVIGATOR_SRC_QUEUE_HPP_
#define NAVIGATOR_SRC_QUEUE_HPP_

template <class T>

class Queue {
 private:
  class Point {
   public:
    T value;
    Point* prev;

    inline Point(T value) : value(value) {}
  };

  Point* first;
  Point* last;
  size_t size;

 public:
  inline size_t getSize() { return size; }

  inline void init() { size = 0; }  // инициализация

  inline void push(const T& value) {  // заполнение
    Point* point = new Point(value);

    switch (size) {  // создание нулевого
      case 0:
        first = point;
        break;

      case 1:  // создание первого
        first->prev = point;
        last = point;
        break;

      default:  // создание второго и всех остальных
        last->prev = point;
        last = point;
        break;
    }

    size += 1;
  }

  inline T pop() {
    T value = first->value;

    Point* point = first;
    first = first->prev;
    size -= 1;

    delete point;
    return value;
  }

  inline const T& peek() const {
    return first->value;
  }  //  значение, что в первом эл-те очереди лежит

  inline bool isEmpty() const { return size == 0; }  // проверка, что не пустой
};
#endif  // NAVIGATOR_SRC_QUEUE_HPP_