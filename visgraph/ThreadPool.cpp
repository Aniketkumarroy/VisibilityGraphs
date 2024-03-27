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
        void start(const char* s){
            m_name = s;
            m_StartPoint = std::chrono::high_resolution_clock::now();
        }

        Timer(const char* s){
            start(s);
        }
        
        Timer(){
            start("");
        }

        void stop(){
            m_EndPoint = std::chrono::high_resolution_clock::now();
            auto start = std::chrono::time_point_cast<std::chrono::microseconds>(m_StartPoint).time_since_epoch().count();
            auto end = std::chrono::time_point_cast<std::chrono::microseconds>(m_EndPoint).time_since_epoch().count();
            std::cout<<m_name<<"\nExecution Time: "<<(end - start)*0.001<<" ms("<<(end - start)<<"us)\n";
        }

        ~Timer(){
            stop();
        }
};

class ThreadPool
{
    public:
        using Task = std::function<void()>;
    
        explicit ThreadPool(std::size_t numThreads){
            start(numThreads);
        }
    
        ~ThreadPool(){
            stop();
        }
    
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
    
        void start(std::size_t numThreads){
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
    
        void stop() noexcept{
            {
                std::unique_lock<std::mutex> lock{mEventMutex};
                mStopping = true;
            }
    
            mEventVar.notify_all();
    
            for (auto &thread : mThreads)
                thread.join();
        }
};

#endif THREADPOOL