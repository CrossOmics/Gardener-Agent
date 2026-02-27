
import asyncio
from concurrent.futures import ThreadPoolExecutor
import multiprocessing
from loguru import logger

class TaskExecutor:
    def __init__(self, max_workers: int = 5):
        # Determine the number of workers based on CPU cores, but not more than max_workers
        cpu_cores = multiprocessing.cpu_count()
        self.num_workers = min(cpu_cores, max_workers)
        logger.info(f"Initializing thread pool with {self.num_workers} workers (CPU cores: {cpu_cores}, Max limit: {max_workers})")
        self.executor = ThreadPoolExecutor(max_workers=self.num_workers)

    async def run_in_thread(self, func, *args, **kwargs):
        """
        Runs a synchronous function in a separate thread to avoid blocking the asyncio event loop.
        """
        loop = asyncio.get_running_loop()
        return await loop.run_in_executor(self.executor, lambda: func(*args, **kwargs))

    def shutdown(self):
        """
        Shuts down the thread pool executor.
        """
        logger.info("Shutting down thread pool executor...")
        self.executor.shutdown(wait=True)

# Global instance
task_executor = TaskExecutor()
