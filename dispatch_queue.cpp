//
// Created by David Winderl on 5/10/21.
//

#include "dispatch_queue.h"


void dispatch_queue_t::thread_handler() {
    while (true) {
        std::unique_lock<std::mutex> lock(_pending_tasks_mutex);

        // wait unless exiting or there's pending tasks
        _task_pending.wait(lock, [this]() {
            return _done || !_pending_tasks.empty();
        });

        if (_done) {
            return;
        }


        try {
            auto task = std::move(_pending_tasks.front());
            _pending_tasks.pop_front();
            lock.unlock();

            task();
        } catch (...) {
            // TODO: I'm not sure what to do here
        }

        _task_completion.notify_one();
    }
}


void dispatch_queue_t::cancel() {
    std::unique_lock<std::mutex> lock(_pending_tasks_mutex);
    _pending_tasks.clear();
    lock.unlock();
}


void dispatch_queue_t::wait() {
    std::unique_lock<std::mutex> lock(_pending_tasks_mutex);
    // block until no pending tasks remain and all workers are idle
    _task_completion.wait(lock, [this]() {
        return _pending_tasks.empty() && _num_busy == 0;
    });
}