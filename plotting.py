import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

# 3D with Constant N
df = pd.read_csv('C:\\Users\\ayush\\OneDrive\\Desktop\\_\\Codes\\CP 4 Lab\\C++\\Fermi-Pasta-Ulam-Tsingou\\FPUT_0.100000.csv')
# use "FPUT_1.000000.csv" for A = 0.1
plt.plot(df['C_time'], df['ET1'], label = "Energy of mode 1")
plt.plot(df['C_time'], df['ET2'], label = "Energy of mode 2")
plt.plot(df['C_time'], df['ET3'], label = "Energy of mode 3")
plt.xlabel('Time')
plt.ylabel('Energy')
plt.title('Energy Evolution of First Three Modes')
plt.legend()
plt.show()

N = 32
# Extract normalized energy columns
norm_columns = [f'norm_E{i}' for i in range(1, N+1)]
normalized_data = df[norm_columns].values

# Create heatmap
plt.figure(figsize=(15/1.5, 8/1.5))
plt.imshow(normalized_data.T, aspect='auto', origin='lower',
           extent=[df['C_time'].iloc[0], df['C_time'].iloc[-1], 1, N],
           cmap='hot', vmin=0, vmax=1)

plt.colorbar(label='Normalized Energy')
plt.xlabel('Time')
plt.ylabel('Mode Index k')
plt.title(f'Normalized Mode Energy Distribution')
plt.show()

# Create meshgrid for surface plot
T, K = np.meshgrid(df['C_time'], range(1, 33))  # Time and mode indices
Z = normalized_data.T  # Energy values (modes x time)

fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d')
surf = ax.plot_surface(T, K, Z, cmap='hot', alpha=1)
ax.set_xlabel('Time')
ax.set_ylabel('Mode k')
ax.set_zlabel('Energy')
ax.set_title('FPUT: 3D Energy Distribution')
fig.colorbar(surf, label='Normalized Energy')  # Only use this line for colorbar
plt.show()

N = 32
# Extract normalized energy columns
RTK_coloumn = [f'RK_t{i}' for i in range(1, N+1)]
normalized_data = df[RTK_coloumn].values

# Create heatmap
plt.figure(figsize=(15/1.5, 8/1.5))
plt.imshow(normalized_data.T, aspect='auto', origin='lower',
           extent=[df['C_time'].iloc[0], df['C_time'].iloc[-1], 1, N],
           cmap='inferno', vmin=0, vmax=1)
plt.xlim(0, 500)
plt.colorbar(label='Recurance')
plt.xlabel('Time')
plt.ylabel('Mode Index k')
plt.title(f'R(t) distribution')
plt.show()


# Create meshgrid for surface plot
T, K = np.meshgrid(df['C_time'], range(1, 33))  # Time and mode indices
Z = normalized_data.T  # Energy values (modes x time)

fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d')
surf = ax.plot_surface(T, K, Z, cmap='viridis', alpha=1)
ax.set_xlabel('Time')
ax.set_ylabel('Mode k')
ax.set_zlabel('Recurance(t)')
ax.set_title('FPUT: 3D Recurance Distribution')
fig.colorbar(surf, label='Recurance')  # Only use this line for colorbar
plt.show()


plt.figure(figsize=(15/2, 8/2))
plt.scatter(df['C_time'], df['R_time'],  s=1)
plt.xlim(0, 15000)
plt.xlabel("Time")
plt.ylabel("Recurance")
plt.show()

plt.figure(figsize=(15/2, 8/2))
plt.plot(df['C_time'], df['R_time'])
plt.xlim(0, 15000)

plt.xlabel("Time")
plt.ylabel("Recurance")
plt.show()