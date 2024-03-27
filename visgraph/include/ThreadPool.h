#ifndef THREADPOOL
#define THREADPOOL

#include<iostream>
#include<future>
#include<shared_mutex>
#include<vector>
#include<queue>
#include<chrono>

class Timer{
    private:
        std::chrono::time_point<std::chrono::high_resolution_clock> m_StartPoint, m_EndPoint;
        const char* m_name;
    public:
        void start(const char* s);
        Timer(const char* s);
        
        Timer();

        void stop();

        ~Timer();
};

class ThreadPool
{
    public:
        using Task = std::function<void()>;
    
        explicit ThreadPool(std::size_t numThreads);
    
        ~ThreadPool();
    
        template<class T>
        auto enqueue(T&& task)->std::future<decltype(task())>{
            auto wrapper = std::make_shared<std::packaged_task<decltype(task()) ()>>(std::move(task));

            {
                std::unique_lock<std::mutex> lock{mEventMutex};
                mTasks.emplace([=] {(*wrapper)();});
            }

            mEventVar.notify_one();
            return wrapper->get_future();
        }
    
    private:
        std::vector<std::thread> mThreads;
    
        std::condition_variable mEventVar;
    
        std::mutex mEventMutex;
        bool mStopping = false;
    
        std::queue<Task> mTasks;
    
        void start(std::size_t numThreads);
    
        void stop() noexcept;
};

#endif