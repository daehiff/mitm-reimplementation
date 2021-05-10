//
// Created by David Winderl on 5/10/21.
//

#ifndef MITM_ALGORITHM_DISPATCH_QUEUE_H
#define MITM_ALGORITHM_DISPATCH_QUEUE_H

#include "utils.h"

class dispatch_queue_t {
public:

    typedef std::function<void()> task_t;

    explicit dispatch_queue_t() {
        this->_done = false;
        const std::size_t n = std::thread::hardware_concurrency();
        for (auto i = 0U; i < n; ++i) {
            _worker_threads.emplace_back(&dispatch_queue_t::thread_handler, this);
        }
    }

    ~dispatch_queue_t() {
        _done = true;
        _task_pending.notify_all(); // wake threads to check for `done` condition and exit`
        for (auto &t : _worker_threads) {
            if (t.joinable()) {
                t.join();
            }
        }
    }

    template<typename F>
    void enqueue(F &&task) {
        std::unique_lock<std::mutex> lock(_pending_tasks_mutex);
        _pending_tasks.emplace_back(std::forward<F>(task));
        lock.unlock();

        _task_pending.notify_one();
    }

    void cancel();

    void wait();

private:

    void thread_handler();

private:

    std::vector<std::thread> _worker_threads;
    std::deque<task_t> _pending_tasks;
    std::mutex _pending_tasks_mutex;

    std::condition_variable _task_pending;
    std::condition_variable _task_completion;

    bool _done;
    int _num_busy = 0;
};


#endif //MITM_ALGORITHM_DISPATCH_QUEUE_H
