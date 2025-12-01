import tkinter as tk
from tkinter import messagebox
import subprocess
import os

# -----------------------------
# 配置路径
# -----------------------------
QUESTION_PATH = r"E:\ChatDev-macnet1\data\question.txt"
SOLVER_PATH = r"E:\ChatDev-macnet1\creat solver_json.py"
RUN_PY_PATH = r"E:\ChatDev-macnet1\run.py"

# -----------------------------
# 主逻辑函数
# -----------------------------
def start_answer():
    # 1️⃣ 获取输入内容
    question = text_input.get("1.0", tk.END).strip()

    if not question:
        messagebox.showwarning("提示", "请输入问题内容！")
        return

    try:
        # 2️⃣ 写入 question.txt
        with open(QUESTION_PATH, "w", encoding="utf-8") as f:
            f.write(question)
        print(f"已写入: {QUESTION_PATH}")

        # 3️⃣ 运行 creat solver_json.py
        messagebox.showinfo("运行中", "正在运行 solver_json.py，请稍候...")
        subprocess.run(["python", SOLVER_PATH], check=True)

        # 4️⃣ 执行 run.py 命令
        command = ["python", RUN_PY_PATH, "--task", question, "--name", "RefrigerantSafety"]
        print("执行命令:", " ".join(command))
        subprocess.run(command, check=True)

        messagebox.showinfo("完成", "任务执行成功！")

    except subprocess.CalledProcessError as e:
        messagebox.showerror("错误", f"程序执行失败：\n{e}")
    except Exception as e:
        messagebox.showerror("错误", f"发生异常：\n{e}")

# -----------------------------
# 构建UI界面
# -----------------------------
root = tk.Tk()
root.title("智能问题求解器")
root.geometry("600x400")

# 输入框
tk.Label(root, text="请输入问题内容：", font=("微软雅黑", 12)).pack(pady=10)
text_input = tk.Text(root, height=10, width=60, font=("微软雅黑", 10))
text_input.pack(padx=10, pady=5)

# 按钮
start_button = tk.Button(root, text="开始回答", font=("微软雅黑", 12, "bold"),
                         bg="#4CAF50", fg="white", width=12, command=start_answer)
start_button.pack(pady=20)

# 进入主循环
root.mainloop()
