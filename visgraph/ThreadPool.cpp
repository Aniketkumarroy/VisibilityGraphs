#include<iostream>
#include<future>
#include<shared_mutex>
#include<vector>
#include<queue>
#include<chrono>
#include"ThreadPool.h"

void Timer::start(const char* s){
    m_name = s;
    m_StartPoint = std::chrono::high_resolution_clock::now();
}

Timer::Timer(const char* s){
    start(s);
}

Timer::Timer(){
    start("");
}

void Timer::stop(){
    m_EndPoint = std::chrono::high_resolution_clock::now();
    auto start = std::chrono::time_point_cast<std::chrono::microseconds>(m_StartPoint).time_since_epoch().count();
    auto end = std::chrono::time_point_cast<std::chrono::microseconds>(m_EndPoint).time_since_epoch().count();
    std::cout<<m_name<<"\nExecution Time: "<<(end - start)*0.001<<" ms("<<(end - start)<<"us)\n";
}

Timer::~Timer(){
    stop();
}


ThreadPool::ThreadPool(std::size_t numThreads){
    start(numThreads);
}

ThreadPool::~ThreadPool(){
    stop();
}

void ThreadPool::start(std::size_t numThreads){
    for (auto i = 0u; i < numThreads; ++i){
        mThreads.emplace_back([=] {
            while (true){
                Task task;

                {
                    std::unique_lock<std::mutex> lock{mEventMutex};

                    mEventVar.wait(lock, [=] { return mStopping || !mTasks.empty(); });

                    if (mStopping && mTasks.empty())
                        break;

                    task = std::move(mTasks.front());
                    mTasks.pop();
                }

                task();
            }
        });
    }
}

void ThreadPool::stop() noexcept{
    {
        std::unique_lock<std::mutex> lock{mEventMutex};
        mStopping = true;
    }

    mEventVar.notify_all();

    for (auto &thread : mThreads)
        thread.join();
}